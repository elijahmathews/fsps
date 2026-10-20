# Cutover to libopensps

Tools for moving this repository into `elijahmathews/libopensps` (`doc/ROADMAP.md`, "Cutover to
libopensps"). libopensps already has its own `main` (README with the logo, LICENSE.md, AGENTS.md,
CONTRIBUTING.md, SECURITY.md, CODE_OF_CONDUCT.md, assets/; Conventional Commits, squash merges,
linear history). The fork arrives as one commit on top of it, followed by the rename. The
sequence for the day itself is in [CHECKLIST.md](CHECKLIST.md).

| File | Purpose |
|---|---|
| `make_import_commit.sh <commit> <libopensps-dir>` | The import: one commit `feat: import elijahmathews/fsps at commit <hash>` on a branch `import` from libopensps's `main`; no push |
| `rename_to_opensps.sh <repo-dir>` | The rename pull request: public C surface and build artifacts to `opensps`, version 0.1.0; idempotent, does not commit |
| `authors.mailmap` | Identities merged when `AUTHORS` is generated (with the evidence for each) |
| `templates/CITATION.cff` | Added by the import (cites the FSPS papers and libopensps) |
| `proposals/README.md` | Sections proposed for libopensps's README (not applied: the repository's README wins) |

Both scripts print `--help`.

## The import commit

`make_import_commit.sh` exports the fork's tree with `git archive <commit>` (so nothing untracked
or ignored, such as build directories or converted data roots, can get in), creates the branch
`import` from `origin/main` of a libopensps clone, copies the tree over it and commits once:

```
feat: import elijahmathews/fsps at commit <hash>

Source: https://github.com/elijahmathews/fsps (tag fsps-final), commit <hash>,
the modernized FSPS fork. Its tree is imported without its history, except:
- left out: data/ (1.4GB), the libraries: published as data releases of
  elijahmathews/opensps-data;
- left out: tests/data/sps_ref_*.bin and tests/data/param_ref/ (59MB),
  the regression references: part of the data release;
- not imported, this repository's own files kept instead: AGENTS.md,
  CLAUDE.md, LICENSE.md, README.md;
- CI workflows run on main (they ran on master and develop);
- added: AUTHORS, CITATION.cff.

History and authorship are recorded in AUTHORS, derived from the fork's
history; the full history stays in elijahmathews/fsps.
```

The author and committer are your git identity; the commit is unsigned and carries no trailers
(`--trailer` adds some).

**libopensps's files win.** Every path that exists on `main` keeps its version there: the fork's
`README.md`, `LICENSE.md` and `AGENTS.md` (and `CLAUDE.md`, a symbolic link to `AGENTS.md` on
`main`) are not imported. The fork's LICENSE carries the same FSPS notice as libopensps's; what
its README and AGENTS.md hold that libopensps's lack is in `proposals/README.md` and in the
fork's AGENTS.md (whose build, test and correctness sections libopensps's AGENTS.md marks as
TODO), to be merged by hand. The fork's `GEMINI.md` (`@AGENTS.md`) comes along and points to
libopensps's AGENTS.md.

**Left out:**

| Path | Size | Where it goes |
|---|---|---|
| `data/` | 1.4 GB | the data releases of elijahmathews/opensps-data, fetched by `tools/fetch_data.py` |
| `tests/data/sps_ref_*.bin` (5 files) | 52 MB | data-v2: the library `references/sps` |
| `tests/data/param_ref/` (35 files) | 7.1 MB | data-v2: the library `references/param` |

The references are copied with `SHA256SUMS` to `<libopensps-dir>.references` (`--refs-dir`) so
that they can be compared with the release; `--keep-refs` imports them instead. Why they leave:
they are generated outputs that every intended physics change regenerates, and each regeneration
would add ~60 MB to the repository's history for good. In the data release they are versioned
with the data they were computed from; the regression tests find them through the data root and
fall back to `tests/data/` (PR #30).

What remains above 1 MB (the script lists it on every run): the test library's MILES fluxes
(`tests/data/fsps-test-data/.../miles/flux.bin`, 2.4 MB), which the unit suite runs on.

**CI filters:** the imported workflows ran on `master` and `develop`; the import points them at
`main` (`--ci-branch`), so that CI runs on the import pull request itself.

**AUTHORS** is generated from the history of `<commit>` (`git log --use-mailmap`, ordered by
number of commits): FSPS contributors are those with commits up to the fork point
(`--fsps-base`, the merge base with `upstream/master`, `82a873508d50`), fork contributors those
after it. Identities are merged only where they are certainly one person (`authors.mailmap`:
Charlie Conroy under six addresses; Benjamin Johnson, Nell Byler and Elijah Mathews under two
each). Two committers appear only as addresses; AUTHORS lists them as "contributors from the
Google Code era of FSPS, known only by these handles" (`tomorrowsville`, `avillaum`). Bots and AI
agents are not listed as authors; one sentence notes that the fork was written with AI agents.

## How the import lands

`main` is protected and squash-merged. Two ways:

1. **Recommended: a pull request from the `import` branch.** Squash-merging a pull request of
   one commit keeps it as one commit with its message (check the message in the merge dialog;
   GitHub fills it from the commit for one-commit pull requests). CI runs on the pull request
   (the imported workflows now trigger on `main`), the protection rules stay in force, and the
   import is reviewable like any other change. Set the repository variable `FSPS_DATA_RELEASE`
   first, so that CI fetches the data and references.
2. **A direct push** of the commit to `main` by you, bypassing the protection: the commit
   arrives byte-for-byte as built (author, date, message), but without CI before it lands.

The rename follows as its own pull request: `refactor!: rename the public C API and library to
opensps`, with a `BREAKING CHANGE:` footer (CHECKLIST.md, step 5).

## Rename scope

**Implemented: the public surface; internal Fortran names kept.**

| Renamed | From | To |
|---|---|---|
| Library version and soname | 3.2.0, `libfsps.so.3` | 0.1.0, `libopensps.so.0` (`--version X.Y.Z`; `--keep-version`) |
| Header | `include/fsps.h` | `include/opensps.h` |
| pkg-config | `fsps.pc` | `opensps.pc` |
| C symbols (all exported `bind(C)` names and header functions) | `fsps_*` | `opensps_*` |
| Header macros (status codes, API version, guard, legacy switch) | `FSPS_*` | `OPENSPS_*` |
| Meson | project `fsps`, `library('fsps')`, `fsps_dep`, `override_dependency('fsps')` | `opensps`, `opensps_dep` |
| Packaging | `fsps.spec`, `libfsps3`, `libfsps-dev` | `opensps.spec`, `libopensps0`, `libopensps-dev` |
| Docs and tests | C tests, `doc/FSPS_C_API.md`, mentions in docs | `opensps_*`, `doc/C_API.md` |
| CI | branch filters | `[main]` (`--branch`; already done by the import) |

The rename never edits README.md, LICENSE.md, AGENTS.md, CONTRIBUTING.md, SECURITY.md,
CODE_OF_CONDUCT.md, CITATION.cff, AUTHORS or assets/.

After the rename the library exports no C function starting with `fsps` (checked with `nm -D`).
What remains of `fsps` in its symbol table are Fortran symbols: module procedures
(`__fsps_*_MOD_*`, namespaced by module) and one external procedure,
`fsps_scenario_batch_of_one_` (`src/physics/fsps_scenario_batch.f90`, PR #22), which sits
outside a module and so is exported under its plain name.

Kept, on purpose:

- **Fortran module, procedure, constant and file names (`fsps_*`)**: not an interface (the
  `.mod` files are not installed). See "Renaming the internal names" for what renaming them
  would cost; libopensps's AGENTS.md currently assumes they are renamed
  (`src/core/opensps_precision.f90`).
- **Environment variables** (`SPS_HOME`, `FSPS_DATA_HOME`, `FSPS_TEST_*`, ...) and the data
  directory names (`share/fsps`, `~/.local/share/fsps`): existing setups keep working;
  `OPENSPS_DATA_HOME`, the preferred name, is read first (PR #29).
- **Internal preprocessor switches** (`FSPS_PINNED`, `FSPS_CUBLAS`).

**Compatibility:** `include/fsps.h` becomes a deprecated header that includes `opensps.h` and
`#define`s each old function and macro name to the new one, generated once from the header as
it was before the rename and frozen afterwards. It is installed alongside `opensps.h`, costs one
file, and keeps C code written against the fork compiling. There is no binary compatibility
(the soname changes anyway). Remove it at the next major version.

**Rerunnable:** the script finds the names to rename in the sources (bind names, header
declarations, header macros), not in a fixed list. Run again on a tree that gained `fsps_*`
entry points (a ported patch), it renames just those; on a renamed tree it changes nothing.
`--check-only` lists what still carries an old name.

**Versions:** libopensps starts at 0.1.0, soname `libopensps.so.0`: the inherited 3.2.0 and
`.so.3` were libfsps's, and a new library name starts its own numbering (before 1.0, an
incompatible change raises the minor version and keeps soname 0; from 1.0 on, the soname follows
the major version). The RPM and Debian package versions follow (`libopensps0`). The C API version
(`OPENSPS_API_VERSION_*`, `opensps_get_driver_version`, now 1.2) stays separate: it describes the
interface the bindings check, changes only when entry points do, and would otherwise jump with
every library release. The `FSPS Version: 3.2` line of the verbose banner also stays: it names
the FSPS version whose physics the code follows.

## Renaming the internal names (open question)

Renaming the internal names too (`src/**/fsps_*.f90` -> `opensps_*.f90`, modules `fsps_*` ->
`opensps_*`, procedures and constants with the prefix, test files `test_fsps_*`) would touch, as
of PR #30: 67 file names and about 2,400 identifier occurrences in 113 files (of ~73,000 lines in
`src/` and `tests/`). Costs:

- **Porting fork patches:** `git am` of a fork patch onto the renamed tree fails (paths and
  identifiers differ). The replay procedure below still works, since the rename is scripted:
  apply on the import commit, rerun the rename stages, carry the difference over.
- **Blame:** within libopensps every line touched would point at the rename commit unless that
  commit is listed in `.git-blame-ignore-revs` (which GitHub honors). The history before the
  import is not in libopensps anyway, so nothing is lost across the cutover beyond that.
- **Diff size:** a large but purely mechanical pull request; review is checking that the tests
  and the regression outputs are unchanged.
- **Pitfalls:** names that must not change (environment variables `FSPS_*`, the data paths
  `tests/data/fsps-test-data` and `tools/fsps_data`, everything under `tools/legacy_parity`,
  which refers to legacy FSPS) and code that derives names from module names (the unit test
  names in `meson.build` cut the `fsps_` prefix with `substring(5)`).

`rename_to_opensps.sh` could do it as an optional second stage (`--internal`): `git mv` of the
files, then `\bfsps_(\w+)` -> `opensps_\1` in Fortran sources, `meson.build` and docs, with the
exclusions above, verified by the unit and regression suites bit for bit. Not implemented;
decide whether the internal names should follow.

## Porting a patch from this repository after the rename

The import commit holds the fork's tree apart from the left-out files, so a fork patch applies
on top of it and the rename can be replayed:

```sh
cd libopensps
imp=$(git log --format=%H --grep='^feat: import elijahmathews/fsps' -1 main)
git switch -c port-A "$imp" && git am /path/to/0001-*.patch
tools/cutover/rename_to_opensps.sh . && git add -A && git -c commit.gpgsign=false commit -m "chore: rename (port)"
git switch -c rename-only "$imp" && tools/cutover/rename_to_opensps.sh . && git add -A && git -c commit.gpgsign=false commit -m "chore: rename"
git switch -c port/A main && git diff rename-only port-A | git apply -3
```
