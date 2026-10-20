#!/usr/bin/env bash
# Rename the public surface of the library from fsps to opensps, in a libopensps checkout.
#
#   tools/cutover/rename_to_opensps.sh [options] <repo-dir>
#
# Rerunnable and idempotent: a second run changes nothing, and a run on a tree that gained new
# fsps_* entry points (e.g. a patch ported from elijahmathews/fsps) renames only those. It edits
# the working tree (and moves files with git mv when <repo-dir> is a git repository); it does not
# commit. What it renames, and what it deliberately keeps (tools/cutover/README.md has the why):
#
#   renamed  library libfsps -> libopensps (soname), header fsps.h -> opensps.h, pkg-config
#            fsps -> opensps, C symbols fsps_* -> opensps_* (every bind(C) name and every
#            function of the public header), public header macros FSPS_* -> OPENSPS_*, Meson
#            project / library / dependency names, packaging (RPM, Debian), C tests, the C API
#            document, CI branch filters, README / CITATION.cff / CONTRIBUTING.md (templates/)
#   kept     Fortran module, procedure, constant and file names (fsps_*), environment variables
#            (SPS_HOME, FSPS_DATA_HOME, FSPS_TEST_*, ...), the installed data directory
#            (share/fsps), internal preprocessor switches (FSPS_PINNED, FSPS_CUBLAS)
#   added    include/fsps.h, a deprecated header that includes opensps.h and maps the old names
#            (written once, from the names present at the first run; frozen afterwards)
#
# Options:
#   --version X.Y.Z     also set the Meson project version (and soversion X); default: unchanged
#   --branch NAME       default branch for CI triggers and AGENTS.md (default: main)
#   --no-templates      do not install README.md / CITATION.cff / CONTRIBUTING.md from templates/
#   --check-only        change nothing; list what still carries the old public names (exit 1 if any)
#   -h, --help          this text
set -euo pipefail

die() { echo "rename_to_opensps: $*" >&2; exit 1; }
usage() { sed -n '2,/^set -euo/p' "$0" | sed -e '/^set -euo/d' -e 's/^# \{0,1\}//'; }

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
version=""
branch="main"
templates=1
check_only=0
args=()
while [ $# -gt 0 ]; do
    case "$1" in
        --version) version=$2; shift 2 ;;
        --branch) branch=$2; shift 2 ;;
        --no-templates) templates=0; shift ;;
        --check-only) check_only=1; shift ;;
        -h|--help) usage; exit 0 ;;
        -*) die "unknown option $1 (see --help)" ;;
        *) args+=("$1"); shift ;;
    esac
done
[ ${#args[@]} -eq 1 ] || { usage >&2; exit 2; }
root=$(cd "${args[0]}" && pwd)
cd "$root"
[ -f meson.build ] || die "$root has no meson.build"
grep -qE "^project\('(fsps|opensps)'" meson.build || die "$root is not the fsps / libopensps source tree"
[ -z "$version" ] || [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || die "--version must be X.Y.Z"

is_git=0
git rev-parse --is-inside-work-tree >/dev/null 2>&1 && is_git=1

# Never edited: this tool and its templates, the credit and license files (they name
# elijahmathews/fsps on purpose), the legacy-FSPS comparison tools (their "fsps" is legacy FSPS),
# the compatibility header, and data.
excluded_re='^(tools/cutover/|tools/legacy_parity/|tools/benchmark_legacy/|AUTHORS$|LICENSE\.md$|include/fsps\.h$|data/|tests/data/)'
compat_marker='libopensps compatibility header'
tmp_list=$(mktemp)
trap 'rm -f "$tmp_list"' EXIT

list_files() {
    if [ $is_git -eq 1 ]; then
        git ls-files
    else
        find . -type f -not -path './.git/*' -printf '%P\n'
    fi | grep -v -E "$excluded_re" | while IFS= read -r f; do
        [ -f "$f" ] && grep -Iq . "$f" 2>/dev/null || continue
        # README / CITATION.cff / CONTRIBUTING.md already written for libopensps are left alone
        case "$f" in README.md|CITATION.cff|CONTRIBUTING.md) grep -q libopensps "$f" && continue ;; esac
        printf '%s\n' "$f"
    done
    return 0
}

move() {  # move $1 to $2 (git mv in a repository); no-op once moved
    [ -e "$1" ] || return 0
    [ -e "$2" ] && die "both $1 and $2 exist"
    mkdir -p "$(dirname "$2")"
    if [ $is_git -eq 1 ] && git ls-files --error-unmatch "$1" >/dev/null 2>&1; then
        git mv "$1" "$2"
    else
        mv "$1" "$2"
    fi
    echo "  moved $1 -> $2"
}

# C symbols still named fsps_*: bind(C) names in the Fortran sources, functions declared in the
# public headers, non-static functions defined in the library's C sources
discover_symbols() {
    {
        find src -name '*.f90' -o -name '*.F90' | xargs -r perl -ne \
            'print "$1\n" while /bind\s*\(\s*c\s*,\s*name\s*=\s*["\x27](fsps_\w+)["\x27]/gi'
        for h in include/*.h; do
            [ -f "$h" ] || continue
            grep -q "$compat_marker" "$h" || perl -ne 'print "$1\n" while /\b(fsps_\w+)\s*\(/g' "$h"
        done
        find src -name '*.c' | xargs -r perl -ne \
            'print "$1\n" if /^(?!\s*static\b)[A-Za-z_][\w\s\*]*?\b(fsps_\w+)\s*\(/'
    } | sort -u
}

# Old public names left anywhere outside the exclusions (exit status 1 if any)
check() {
    local hits left
    hits=$(list_files | xargs -r -d '\n' grep -nIE \
        "#include *[<\"]fsps\.h[>\"]|\blibfsps|\bpkgconfig\(fsps\)|\bfsps\.pc\b|^project\('fsps'|library\('fsps'|override_dependency\('fsps'|dependency\('fsps'\)|\bfsps_dep\b" \
        2>/dev/null || true)
    left=$(discover_symbols)
    [ -z "$left" ] || hits+=$'\n'"C symbols still named fsps_*: $(printf '%s ' $left)"
    if [ $is_git -eq 1 ]; then
        local paths
        paths=$(git ls-files | grep -E '(^|/)(libfsps[^/]*|fsps\.spec|FSPS_C_API\.md)$' || true)
        [ -z "$paths" ] || hits+=$'\n'"Old file names: $(printf '%s ' $paths)"
    fi
    hits=$(printf '%s\n' "$hits" | grep . || true)
    if [ -n "$hits" ]; then
        echo "Old public names remain:"
        printf '%s\n' "$hits" | sed 's/^/  /'
        return 1
    fi
    return 0
}

if [ $check_only -eq 1 ]; then
    if check; then echo "No old public names left."; exit 0; else exit 1; fi
fi

public_header=include/opensps.h
if [ -f include/fsps.h ] && ! grep -q "$compat_marker" include/fsps.h; then
    public_header=include/fsps.h
fi
[ -f "$public_header" ] || die "no public header (include/fsps.h or include/opensps.h)"

symbols=$(discover_symbols)
# Macros of the public header (status codes, API version, include guard, legacy switch)
macros=$(perl -ne 'print "$1\n" if /^\s*#\s*(?:define|ifdef|ifndef|endif\s*\/\/|if\s+defined)\s*\(?\s*(FSPS_\w+)/' \
    "$public_header" | sort -u)

# A Fortran procedure bound to C without name= would keep an fsps_ symbol: refuse
unnamed=$(grep -rniE '^[^!]*\b(subroutine|function)\s+fsps_\w+.*bind\s*\(\s*c\s*\)' src || true)
[ -z "$unnamed" ] || die "bind(C) without name= on fsps_ procedures (give them a name):"$'\n'"$unnamed"

echo "== Renaming in $root"
echo "  $(printf '%s\n' "$symbols" | grep -c . || true) C symbols and $(printf '%s\n' "$macros" | grep -c . || true) public header macros named fsps_* / FSPS_*"

# ---------------------------------------------------------------------------------------------
# 1. Files that move
# ---------------------------------------------------------------------------------------------
old_header_text=""
if [ "$public_header" = include/fsps.h ]; then
    old_header_text=$(cat include/fsps.h)
    move include/fsps.h include/opensps.h
fi
move doc/FSPS_C_API.md doc/C_API.md
move packaging/fsps.spec packaging/opensps.spec
for f in debian/libfsps*.install; do
    [ -e "$f" ] && move "$f" "${f/libfsps/libopensps}"
done

list_files > "$tmp_list"

# ---------------------------------------------------------------------------------------------
# 2. C symbols and public header macros (exact identifiers). In Fortran only inside bind(C)
#    names and on preprocessor lines: Fortran procedure and constant names stay.
# ---------------------------------------------------------------------------------------------
export SYM_RE=$(printf '%s\n' "$symbols" | grep . | paste -sd'|' || true)
export MAC_RE=$(printf '%s\n' "$macros" | grep . | paste -sd'|' || true)
while IFS= read -r f; do
    case "$f" in
        *.f90|*.F90)
            perl -i -pe '
                my $m = $ENV{MAC_RE};
                s/(name\s*=\s*["\x27])fsps_(\w+["\x27])/${1}opensps_$2/gi if /bind\s*\(\s*c/i;
                if ($m ne "" && /^\s*#/) { s/\b(?:$m)\b/OPENSPS_$&/g; s/OPENSPS_FSPS_/OPENSPS_/g; }
            ' "$f" ;;
        *)
            # (not in file names such as src/core/fsps_setup_lock.c, nor on Markdown lines about
            # Fortran imports, where the same name is a Fortran procedure)
            MD=$([[ "$f" == *.md ]] && echo 1 || echo 0) perl -i -pe '
                my $s = $ENV{SYM_RE}; my $m = $ENV{MAC_RE};
                if ($s ne "" && !($ENV{MD} && /\bimports?\b/)) {
                    s/(?<![\w\/.])\b(?:$s)\b(?!\.(?:c|h|f90|F90|o|mod)\b)/opensps_$&/g;
                    s/opensps_fsps_/opensps_/g;
                }
                if ($m ne "") { s/\b(?:$m)\b/OPENSPS_$&/g; s/OPENSPS_FSPS_/OPENSPS_/g; }
            ' "$f" ;;
    esac
done < "$tmp_list"

# ---------------------------------------------------------------------------------------------
# 3. Header, library, pkg-config and package names in text. fsps.h -> opensps.h except on lines
#    about the deprecated compatibility header and on the packaging lines that install it.
# ---------------------------------------------------------------------------------------------
while IFS= read -r f; do
    perl -0 -i -pe '
        my $has = /opensps\.h/;
        s{^(.*)$}{
            my $l = $1;
            ($l =~ /deprecated/i || ($has && $l =~ m{^(?:%\{_includedir\}/|usr/include/)fsps\.h$}))
                ? $l : ($l =~ s/\bfsps\.h\b/opensps.h/gr)
        }gme;
        s/\blibfsps/libopensps/g;
        s/\bpkgconfig\(fsps\)/pkgconfig(opensps)/g;
        s/\bfsps\.pc\b/opensps.pc/g;
        s/\bFSPS_C_API\.md\b/C_API.md/g;
        s/\bfsps\.spec\b/opensps.spec/g;
        s/\bdependency\(\x27fsps\x27\)/dependency(\x27opensps\x27)/g;
        s/\bfsps_dep\b/opensps_dep/g;
    ' "$f"
done < "$tmp_list"

# ---------------------------------------------------------------------------------------------
# 4. Meson
# ---------------------------------------------------------------------------------------------
perl -0 -i -pe '
    s/^project\(\x27fsps\x27/project(\x27opensps\x27/m;
    s/library\(\x27fsps\x27/library(\x27opensps\x27/g;
    s/\bfsps_dep\b/opensps_dep/g;
    s/override_dependency\(\x27fsps\x27/override_dependency(\x27opensps\x27/g;
    s/install_headers\(\x27include\/opensps\.h\x27\)\n/install_headers(\x27include\/opensps.h\x27, \x27include\/fsps.h\x27)  # fsps.h: deprecated compatibility header\n/;
    s/(pkg\.generate\(libopensps,\s*\n\s*name\s*:\s*)\x27fsps\x27/${1}\x27opensps\x27/;
    s/(pkg\.generate\(libopensps,.*?description\s*:\s*)\x27Flexible Stellar Population Synthesis\x27/${1}\x27Stellar population synthesis library derived from FSPS\x27/s;
    s{(pkg\.generate\(libopensps,.*?url\s*:\s*)\x27https://github.com/cconroy20/fsps\x27}{${1}\x27https://github.com/elijahmathews/libopensps\x27}s;
' meson.build
# A tree without data/ must still install: install the data directory only where it exists (the
# libraries are otherwise fetched into FSPS_DATA_HOME)
if grep -q "^ *install_subdir('data'" meson.build && ! grep -q "is_dir('data')" meson.build; then
    perl -0 -i -pe 's/^( *)(install_subdir\(\x27data\x27[^\n]*\n)/${1}if import(\x27fs\x27).is_dir(\x27data\x27)\n${1}    ${2}${1}endif\n/m' meson.build
fi
if [ -n "$version" ]; then
    perl -0 -i -pe "s/(^project\(\x27opensps\x27.*?version\s*:\s*)\x27[^\x27]*\x27/\${1}\x27$version\x27/ms;
                    s/(soversion\s*:\s*)\x27[^\x27]*\x27/\${1}\x27${version%%.*}\x27/" meson.build
fi
soversion=$(perl -ne 'print $1 if /soversion\s*:\s*\x27([^\x27]*)\x27/' meson.build)
proj_version=$(perl -0 -ne 'print $1 if /^project\(.*?version\s*:\s*\x27([^\x27]*)\x27/ms' meson.build)

# ---------------------------------------------------------------------------------------------
# 5. Packaging (RPM, Debian, the packaging workflow)
# ---------------------------------------------------------------------------------------------
if [ -f packaging/opensps.spec ]; then
    perl -i -pe '
        s/^Name:(\s+)fsps\s*$/Name:${1}opensps\n/;
        s/^Summary:(\s+)Flexible Stellar Population Synthesis \(FSPS\)\s*$/Summary:${1}Stellar population synthesis library derived from FSPS\n/;
        s{^URL:(\s+)https://github.com/elijahmathews/fsps\s*$}{URL:${1}https://github.com/elijahmathews/libopensps\n};
        s/^Flexible Stellar Population Synthesis \(FSPS\) is a stellar population$/libopensps is a stellar population synthesis library derived from FSPS/;
        s/^synthesis code\.$/(Flexible Stellar Population Synthesis, Conroy et al.)./;
        s/# The base \x27fsps\x27 package/# The base \x27opensps\x27 package/;
    ' packaging/opensps.spec
    grep -q '^%{_includedir}/fsps.h$' packaging/opensps.spec ||
        perl -i -pe 's{^(%\{_includedir\}/opensps\.h)\n}{$1\n%{_includedir}/fsps.h\n}' packaging/opensps.spec
    perl -i -pe "s{(%\\{_libdir\\}/libopensps\\.so\\.)[0-9]+\\*}{\${1}${soversion}*}" packaging/opensps.spec
    [ -z "$version" ] || perl -i -pe "s/^Version:(\\s+)\\S+/Version:\${1}$version/" packaging/opensps.spec
fi
if [ -d debian ]; then
    perl -i -pe 's/^Source: fsps$/Source: opensps/;
                 s{^Homepage: https://github.com/elijahmathews/fsps$}{Homepage: https://github.com/elijahmathews/libopensps};
                 s/^Maintainer: FSPS Packager/Maintainer: libopensps Packager/;
                 s/^Description: FSPS shared library$/Description: libopensps shared library/;
                 s/^ Shared library for FSPS \(Flexible Stellar Population Synthesis\)\.$/ Shared library of libopensps, derived from FSPS (Flexible Stellar Population Synthesis)./' debian/control
    perl -i -pe 's/^fsps \(/opensps (/' debian/changelog
    perl -i -pe 's/^Upstream-Name: fsps$/Upstream-Name: opensps/;
                 s{^Source: https://github.com/elijahmathews/fsps$}{Source: https://github.com/elijahmathews/libopensps}' debian/copyright
    perl -i -pe 's{/fsps-\$1}{/opensps-\$1}; s{github.com/elijahmathews/fsps/tags}{github.com/elijahmathews/libopensps/tags}' debian/watch
    # The runtime package is named after the soname: libopensps<soversion>
    for f in debian/libopensps[0-9]*.install; do
        [ -e "$f" ] || continue
        want="debian/libopensps${soversion}.install"
        [ "$f" = "$want" ] || move "$f" "$want"
        perl -i -pe "s{(libopensps\\.so\\.)[0-9]+\\*}{\${1}${soversion}*}" "$want"
    done
    perl -i -pe "s/\\blibopensps[0-9]+\\b/libopensps${soversion}/g" debian/control
    if [ -f debian/libopensps-dev.install ] && ! grep -q '^usr/include/fsps.h$' debian/libopensps-dev.install; then
        perl -i -pe 's{^(usr/include/opensps\.h)\n}{$1\nusr/include/fsps.h\n}' debian/libopensps-dev.install
    fi
fi
if [ -f .github/workflows/packaging.yml ]; then
    perl -i -pe 's/\bfsps-\$VERSION/opensps-\$VERSION/g; s{\.\./fsps_\*}{../opensps_*}g;' .github/workflows/packaging.yml
    perl -i -pe "s/\\blibopensps[0-9]+_\\*/libopensps${soversion}_*/g" .github/workflows/packaging.yml
fi

# ---------------------------------------------------------------------------------------------
# 6. CI triggers and the agent guide: the default branch of libopensps
# ---------------------------------------------------------------------------------------------
for wf in .github/workflows/*.yml; do
    [ -f "$wf" ] || continue
    perl -i -pe "s/^(\\s*branches:\\s*)\\[(?:master, develop|master|master, $branch|$branch, develop)\\]/\${1}[$branch]/" "$wf"
done
if [ -f AGENTS.md ]; then
    perl -i -pe "s/topic branch from \`develop\`/topic branch from \`$branch\`/;
                 s/against \`develop\`/against \`$branch\`/;
                 s/^A modernized fork of FSPS \\(Flexible Stellar Population Synthesis; Conroy et al.\\), in Fortran\$/libopensps: a modernized FSPS (Flexible Stellar Population Synthesis; Conroy et al.), in Fortran/" AGENTS.md
fi

# ---------------------------------------------------------------------------------------------
# 7. The C API document's title, and the compatibility header
# ---------------------------------------------------------------------------------------------
if [ -f doc/C_API.md ]; then
    perl -i -pe 's/^# FSPS C Driver API$/# libopensps C API/;
        s/^This document describes the C driver API exposed by FSPS\./This document describes the C API of libopensps: header `opensps.h`, library `libopensps`, pkg-config `opensps`, prefix `opensps_` (macros `OPENSPS_`)./' doc/C_API.md
    grep -q 'deprecated compatibility header' doc/C_API.md ||
        perl -0 -i -pe 's/(\n## Architecture\n)/\n`fsps.h` (deprecated compatibility header) maps the names of the elijahmathews\/fsps C API (prefix `fsps_`) to these; it is frozen at the rename and will be removed in a future major version.\n$1/' doc/C_API.md
fi

if [ -n "$old_header_text" ] && [ ! -f include/fsps.h ]; then
    old_syms=$(printf '%s\n' "$old_header_text" | perl -ne 'print "$1\n" while /\b(fsps_\w+)\s*\(/g' | sort -u)
    old_macs=$(printf '%s\n' "$old_header_text" |
        perl -ne 'print "$1\n" if /^\s*#\s*define\s+(FSPS_\w+)/' | grep -v '^FSPS_H$' | sort -u)
    {
        echo "/*"
        echo " * fsps.h: $compat_marker (deprecated)."
        echo " *"
        echo " * The C API of libopensps is declared in opensps.h, with the prefixes opensps_ (functions)"
        echo " * and OPENSPS_ (macros). This header maps the names of the elijahmathews/fsps C API to"
        echo " * them, so that existing C sources compile unchanged. It is frozen: entry points added"
        echo " * after the rename have only their opensps_ names. It will be removed in a future major"
        echo " * version; include opensps.h in new code."
        echo " */"
        echo "#ifndef FSPS_COMPAT_H"
        echo "#define FSPS_COMPAT_H"
        echo
        echo "/* The legacy declarations were enabled by FSPS_ENABLE_LEGACY */"
        echo "#if defined(FSPS_ENABLE_LEGACY) && !defined(OPENSPS_ENABLE_LEGACY)"
        echo "#define OPENSPS_ENABLE_LEGACY"
        echo "#endif"
        echo
        echo "#include \"opensps.h\""
        echo
        echo "/* Functions */"
        for s in $old_syms; do printf '#define %-40s opensps_%s\n' "$s" "${s#fsps_}"; done
        echo
        echo "/* Status codes and API version */"
        for m in $old_macs; do printf '#define %-40s OPENSPS_%s\n' "$m" "${m#FSPS_}"; done
        echo
        echo "#endif /* FSPS_COMPAT_H */"
    } > include/fsps.h
    [ $is_git -eq 1 ] && git add include/fsps.h
    echo "  wrote include/fsps.h (deprecated compatibility header: $(printf '%s\n' "$old_syms" | grep -c .) functions, $(printf '%s\n' "$old_macs" | grep -c .) macros)"
fi

# ---------------------------------------------------------------------------------------------
# 8. README.md, CITATION.cff, CONTRIBUTING.md from the templates (unless already libopensps's)
# ---------------------------------------------------------------------------------------------
if [ $templates -eq 1 ] && [ -d "$script_dir/templates" ]; then
    for t in README.md CITATION.cff CONTRIBUTING.md; do
        [ -f "$script_dir/templates/$t" ] || continue
        if [ ! -f "$t" ] || ! grep -q 'libopensps' "$t"; then
            cp "$script_dir/templates/$t" "$t"
            [ $is_git -eq 1 ] && git add "$t"
            echo "  installed $t from tools/cutover/templates"
        fi
    done
fi

echo "== Project $(perl -ne 'print $1 if /^project\(\x27(\w+)\x27/' meson.build) $proj_version, soname libopensps.so.$soversion"
if check; then
    echo "No old public names left. Review with: git -C $root status; git -C $root diff --stat"
else
    exit 1
fi
