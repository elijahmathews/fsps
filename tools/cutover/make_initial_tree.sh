#!/usr/bin/env bash
# Build the first commit of libopensps from one commit of this repository (elijahmathews/fsps).
#
#   tools/cutover/make_initial_tree.sh [options] <fork-commit> <out-dir>
#
# Exports the tree at <fork-commit> (git archive), leaves out data/ (the libraries move to
# versioned data archives fetched by tools/fetch_data.py), writes AUTHORS (the history that
# recorded authorship is not carried over) and adds one line to LICENSE.md for the fork's
# contributors, then initializes a git repository in <out-dir> with a single unsigned commit:
#
#   libopensps: the state of elijahmathews/fsps at commit <full hash>
#
# It never adds a remote and never pushes. See tools/cutover/README.md and CHECKLIST.md.
#
# Options:
#   --repo DIR            the fork's clone (default: the repository holding this script)
#   --author "N <e>"      author and committer of the snapshot (default: git config user.* of --repo)
#   --branch NAME         initial branch (default: main)
#   --fsps-base REV       last commit of upstream FSPS in the fork's history; AUTHORS credits
#                         commits up to it to FSPS and later ones to the fork (default: merge-base
#                         with upstream/master if that remote exists, else 82a873508d50, the fork point)
#   --exclude-refs DIR    also leave out the regression references tests/data/sps_ref_*.bin
#                         (54 MB) and copy them, with SHA256SUMS, to DIR (e.g. for a release asset).
#                         Off by default: see "What is left out" in tools/cutover/README.md.
#   --trailer TEXT        append a trailer line to the commit message (repeatable). None by default:
#                         the commit is a snapshot of other commits, not new work.
#   --date DATE           author/committer date of the snapshot (default: now)
#   -h, --help            this text
set -euo pipefail

die() { echo "make_initial_tree: $*" >&2; exit 1; }
usage() { sed -n '2,/^set -euo/p' "$0" | sed -e '/^set -euo/d' -e 's/^# \{0,1\}//'; }

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=""
author=""
branch="main"
fsps_base=""
refs_dir=""
date=""
trailers=()
default_fsps_base=82a873508d500ca353bbb922459bf928498f7a72

args=()
while [ $# -gt 0 ]; do
    case "$1" in
        --repo) repo=$2; shift 2 ;;
        --author) author=$2; shift 2 ;;
        --branch) branch=$2; shift 2 ;;
        --fsps-base) fsps_base=$2; shift 2 ;;
        --exclude-refs) refs_dir=$2; shift 2 ;;
        --trailer) trailers+=("$2"); shift 2 ;;
        --date) date=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        --) shift; args+=("$@"); break ;;
        -*) die "unknown option $1 (see --help)" ;;
        *) args+=("$1"); shift ;;
    esac
done
[ ${#args[@]} -eq 2 ] || { usage >&2; exit 2; }
rev=${args[0]}
out=${args[1]}

[ -n "$repo" ] || repo=$(git -C "$script_dir" rev-parse --show-toplevel)
git -C "$repo" rev-parse --git-dir >/dev/null 2>&1 || die "$repo is not a git repository"
commit=$(git -C "$repo" rev-parse --verify "${rev}^{commit}") || die "no commit $rev in $repo"

if [ -z "$author" ]; then
    name=$(git -C "$repo" config user.name) || die "no user.name in $repo; pass --author"
    email=$(git -C "$repo" config user.email) || die "no user.email in $repo; pass --author"
    author="$name <$email>"
fi
[[ "$author" =~ ^(.+)\ \<([^>]+)\>$ ]] || die "--author must be 'Name <email>'"
author_name=${BASH_REMATCH[1]}
author_email=${BASH_REMATCH[2]}

if [ -z "$fsps_base" ]; then
    if git -C "$repo" rev-parse --verify -q upstream/master >/dev/null; then
        fsps_base=$(git -C "$repo" merge-base "$commit" upstream/master)
    else
        fsps_base=$default_fsps_base
    fi
fi
fsps_base=$(git -C "$repo" rev-parse --verify "${fsps_base}^{commit}") || die "no commit $fsps_base"
git -C "$repo" merge-base --is-ancestor "$fsps_base" "$commit" || die "$fsps_base is not an ancestor of $commit"

if [ -e "$out" ] && [ -n "$(ls -A "$out" 2>/dev/null)" ]; then
    die "$out exists and is not empty"
fi
mkdir -p "$out"
out=$(cd "$out" && pwd)

echo "== Exporting $commit"
git -C "$repo" archive --format=tar "$commit" | tar -x -C "$out"

# ---------------------------------------------------------------------------------------------
# What is left out
# ---------------------------------------------------------------------------------------------
removed=()
if [ -d "$out/data" ]; then
    removed+=("data/ ($(du -sh "$out/data" | cut -f1); the libraries, distributed as data archives)")
    rm -rf "$out/data"
fi
if [ -n "$refs_dir" ]; then
    mkdir -p "$refs_dir"
    refs_dir=$(cd "$refs_dir" && pwd)
    shopt -s nullglob
    refs=("$out"/tests/data/sps_ref_*.bin)
    shopt -u nullglob
    if [ ${#refs[@]} -gt 0 ]; then
        cp "${refs[@]}" "$refs_dir/"
        (cd "$refs_dir" && sha256sum sps_ref_*.bin > SHA256SUMS)
        removed+=("tests/data/sps_ref_*.bin ($(du -ch "${refs[@]}" | tail -1 | cut -f1); copied to $refs_dir)")
        rm -f "${refs[@]}"
    fi
fi

# Anything else above 1 MB is reported (and kept): the list should hold only known files
echo "== Files above 1 MB in the new tree"
big=$(cd "$out" && find . -type f -size +1024k -printf '%s\t%P\n' | sort -rn)
if [ -n "$big" ]; then
    printf '%s\n' "$big" | awk -F'\t' '{printf "  %6.1f MB  %s\n", $1/1048576, $2}'
else
    echo "  (none)"
fi

# ---------------------------------------------------------------------------------------------
# AUTHORS and the license line for the fork's contributors
# ---------------------------------------------------------------------------------------------
mailmap="$script_dir/authors.mailmap"
[ -f "$mailmap" ] || die "missing $mailmap"

# People (mailmap-resolved "Name <email>") by number of commits over a revision range; bots and
# AI agents are left out (noted separately). The years are those of the git history only (FSPS
# predates it), so they are not printed.
people() {
    git -C "$repo" -c mailmap.file="$mailmap" log --use-mailmap --format='%aN <%aE>%x09%ad' \
        --date=format:%Y "$@" |
        grep -v -E '^(Claude|Google Code Exporter|GitHub) <' |
        awk -F'\t' '{n[$1]++; if (!($1 in lo) || $2 < lo[$1]) lo[$1] = $2; if ($2 > hi[$1]) hi[$1] = $2}
                    END {for (p in n) printf "%d\t%s\t%s\t%s\n", n[p], p, lo[p], hi[p]}' |
        sort -t$'\t' -k1,1nr -k2,2
}
# A named person ("Name <email>" with a real name) or a bare handle (name = an address)
format_people() {
    local named="" handles=""
    while IFS=$'\t' read -r _n who lo hi; do
        local nm=${who% <*} years=$lo
        [ "$hi" != "$lo" ] && years="$lo-$hi"
        if [[ "$nm" == *@* ]]; then
            handles+="${handles:+, }${nm%%@*} (commits $years)"
        else
            named+="- $nm"$'\n'
        fi
    done
    printf '%s' "$named"
    [ -n "$handles" ] && printf '\nCommits recorded under an address only, not yet tied to a name: %s.\n' "$handles"
    return 0
}

fsps_people=$(people "$fsps_base" | format_people)
fork_people=$(people "$fsps_base..$commit" | format_people)
n_ai=$(git -C "$repo" log --format='%aE%n%(trailers:key=Co-authored-by,valueonly)' "$fsps_base..$commit" |
       grep -c -i 'noreply@anthropic.com' || true)
short_base=$(git -C "$repo" rev-parse --short=12 "$fsps_base")

cat > "$out/AUTHORS" <<EOF
# Authors

libopensps is derived from FSPS, Flexible Stellar Population Synthesis, the work of Charlie
Conroy and contributors (https://github.com/cconroy20/fsps). The models are described in
Conroy, Gunn & White (2009, ApJ, 699, 486) and Conroy & Gunn (2010, ApJ, 712, 833); please cite
both papers in work that uses this library (see README.md and CITATION.cff).

This repository starts from a snapshot of the modernized fork elijahmathews/fsps at commit
$commit
without its history. The full history, which records who wrote what, stays available in that
repository (tag fsps-final). The lists below are derived from it, ordered by number of commits.

## FSPS (commit history up to $short_base, where the fork began)

$fsps_people

## The modernized fork (elijahmathews/fsps)

$fork_people
EOF
if [ "$n_ai" -gt 0 ]; then
    cat >> "$out/AUTHORS" <<'EOF'

Much of the fork's code was written with AI coding agents (Claude, by Anthropic) under the
direction and review of the fork's maintainer; the commit messages record them as co-authors.
EOF
fi

# LICENSE.md: the FSPS notice unchanged, one line for the fork's contributors after it
fork_years=$(git -C "$repo" log --format=%ad --date=format:%Y "$fsps_base..$commit" | sort -u | sed -n '1p;$p' | paste -sd- | sed 's/^\([0-9]*\)-\1$/\1/')
fork_line="Copyright © ${fork_years} Elijah Mathews & contributors (the modernized fork, libopensps)."
if [ -f "$out/LICENSE.md" ] && ! grep -qF "the modernized fork, libopensps" "$out/LICENSE.md"; then
    grep -q '^Copyright © 2009-[0-9]* Charlie Conroy & contributors\.$' "$out/LICENSE.md" ||
        die "LICENSE.md does not have the expected FSPS copyright line"
    awk -v line="$fork_line" '{print} /^Copyright © 2009-[0-9]* Charlie Conroy & contributors\.$/ {print line}' \
        "$out/LICENSE.md" > "$out/LICENSE.md.new"
    mv "$out/LICENSE.md.new" "$out/LICENSE.md"
fi

# ---------------------------------------------------------------------------------------------
# The repository and its one commit
# ---------------------------------------------------------------------------------------------
echo "== Creating the repository"
git -C "$out" init -q -b "$branch"
git -C "$out" config user.name "$author_name"
git -C "$out" config user.email "$author_email"
git -C "$out" config commit.gpgsign false
git -C "$out" add -A
msg="libopensps: the state of elijahmathews/fsps at commit $commit"
if [ ${#trailers[@]} -gt 0 ]; then
    msg+=$'\n'
    for t in "${trailers[@]}"; do msg+=$'\n'"$t"; done
fi
date_env=()
[ -n "$date" ] && date_env=(GIT_AUTHOR_DATE="$date" GIT_COMMITTER_DATE="$date")
env "${date_env[@]}" git -C "$out" commit -q --no-verify -m "$msg"
git -C "$out" gc -q --aggressive --prune=now

echo
echo "== Done: $out"
git -C "$out" log --format='  %H%n  %an <%ae>, %ad%n  %s' -1
echo "  files: $(git -C "$out" ls-files | wc -l), checked out: $(du -sh --exclude=.git "$out" | cut -f1), .git: $(du -sh "$out/.git" | cut -f1)"
echo "  left out:"
for r in "${removed[@]}"; do echo "    $r"; done
echo "  added: AUTHORS; LICENSE.md: \"$fork_line\""
echo "No remote was added and nothing was pushed."
