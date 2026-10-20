# Cutover to libopensps

Tools for moving this repository to `elijahmathews/libopensps` (`doc/ROADMAP.md`, "Cutover to
libopensps"). The sequence for the day itself is in [CHECKLIST.md](CHECKLIST.md).

| File | Purpose |
|---|---|
| `make_initial_tree.sh <commit> <out-dir>` | The first commit of libopensps: the tree at `<commit>` without `data/`, plus `AUTHORS` and a license line; one unsigned commit, no remote, no push |
| `rename_to_opensps.sh <repo-dir>` | The rename pull request in libopensps: public C surface and build artifacts to `opensps`; idempotent, does not commit |
| `authors.mailmap` | Identities merged when `AUTHORS` is generated (with the evidence for each) |
| `templates/` | `README.md`, `CITATION.cff`, `CONTRIBUTING.md` for libopensps (installed by the rename) |

Both scripts print `--help`.

## What the first commit leaves out

The snapshot is `git archive <commit>`, so nothing untracked or ignored (build directories,
`.mod` files, converted data roots) can get in. Of the tracked files it leaves out:

| Path | Size | Where it goes |
|---|---|---|
| `data/` | 1.4 GB | the data release (elijahmathews/opensps-data), fetched by `tools/fetch_data.py` |
| `tests/data/sps_ref_*.bin` (5 files) | 52 MB | the data release, as a versioned library of regression references |
| `tests/data/param_ref/` (35 files) | 7.1 MB | the same library |

The references are copied with `SHA256SUMS` to `<out-dir>.references` (`--refs-dir`), to check
them against the release; `--keep-refs` keeps them in the tree instead. Why they leave: they are
generated outputs, and every intended physics change regenerates them; in the repository each
regeneration would add ~60 MB to the history for good.
In the data release they are versioned with the data they were computed from, and
the regression tests find them through the data root (falling back to `tests/data/`).

What remains above 1 MB (the script lists it on every run): the test library's MILES fluxes
(`tests/data/fsps-test-data/.../miles/flux.bin`, 2.4 MB), which the unit suite runs on. Also
kept, all small: `doc/MANUAL.pdf` (the FSPS manual), `tests/data/param_cases.psv` (the case
matrix), editor and agent files (`.vscode/`, `.fortls`, `CLAUDE.md`, `GEMINI.md`).

## AUTHORS and the license

The history that records authorship is not carried over, so the snapshot adds `AUTHORS`,
generated from the history of `<commit>` (`git log --use-mailmap`, ordered by number of
commits): FSPS contributors are those with commits up to the fork point (`--fsps-base`, the
merge base with `upstream/master`, `82a873508d50`), fork contributors those after it. Identities
are merged only where they are certainly one person (`authors.mailmap`: Charlie Conroy under six
addresses, Benjamin Johnson, Nell Byler, Elijah Mathews under two each). Two committers appear
only as addresses; `AUTHORS` lists them as "contributors from the Google Code era of FSPS,
known only by these handles" (`tomorrowsville`, `avillaum`), without guessing names. Bots and AI
agents are not listed as authors; a sentence notes that the fork was written with AI agents.

`LICENSE.md` keeps the FSPS notice unchanged and gains one line after it:
`Copyright © 2026 Elijah Mathews & contributors (the modernized fork, libopensps).`

## Rename scope

**Recommendation: rename the public surface, keep the internal Fortran names.** Implemented:

| Renamed | From | To |
|---|---|---|
| Library version and soname | 3.2.0, `libfsps.so.3` | 0.1.0, `libopensps.so.0` (`--version X.Y.Z`; `--keep-version`) |
| Header | `include/fsps.h` | `include/opensps.h` |
| pkg-config | `fsps.pc` | `opensps.pc` |
| C symbols (all exported `bind(C)` names and header functions) | `fsps_*` | `opensps_*` |
| Header macros (status codes, API version, guard, legacy switch) | `FSPS_*` | `OPENSPS_*` |
| Meson | project `fsps`, `library('fsps')`, `fsps_dep`, `override_dependency('fsps')` | `opensps`, `opensps_dep` |
| Packaging | `fsps.spec`, `libfsps3`, `libfsps-dev` | `opensps.spec`, `libopensps3`, `libopensps-dev` |
| Docs and tests | C tests, `doc/FSPS_C_API.md`, mentions in docs | `opensps_*`, `doc/C_API.md` |
| CI | branch filters `[master, develop]` | `[main]` (`--branch`) |

After the rename the library exports no C-linkage symbol starting with `fsps` (checked with
`nm -D`); the only `fsps` in its symbol table are Fortran module symbols (`__fsps_*_MOD_*`).

Kept, on purpose:

- **Fortran module, procedure, constant and file names (`fsps_*`).** They are not part of any
  interface: the `.mod` files are not installed, and the symbols are namespaced by module.
  Renaming them would touch nearly every line that `use`s a module and every file name, which
  breaks `git am` of the patches still to be ported from this repository and makes `git blame`
  across the cutover harder, for no change a user can see. If wanted, it can be a separate,
  purely mechanical pull request once the porting window has closed.
- **Environment variables** (`SPS_HOME`, `FSPS_DATA_HOME`, `FSPS_TEST_*`, ...) and the data
  directory names (`share/fsps`, `~/.local/share/fsps`): existing setups and CI keep working.
  `OPENSPS_DATA_HOME`, the preferred name, is already read first (PR #29).
- **Internal preprocessor switches** (`FSPS_PINNED`, `FSPS_CUBLAS`): build internals.

**Compatibility:** `include/fsps.h` becomes a deprecated header that includes `opensps.h` and
`#define`s each old function and macro name to the new one, generated once from the header as
it was before the rename and frozen afterwards. It is installed alongside `opensps.h`. It costs
one file and keeps any C code written against the fork compiling unchanged; there is no binary
compatibility (the soname changes anyway). Remove it at the next major version.

**Rerunnable:** the script discovers the names to rename from the sources (bind names, header
declarations, header macros) rather than from a fixed list, so running it again on a tree that
gained `fsps_*` entry points (a ported patch) renames just those, and running it on a renamed
tree changes nothing. `--check-only` lists what still carries an old name.

**Versions:** libopensps starts at 0.1.0, soname `libopensps.so.0`: the inherited 3.2.0 and
`.so.3` were libfsps's, and a new library name starts its own numbering (before 1.0, an
incompatible change raises the minor version and keeps soname 0; from 1.0 on, the soname
follows the major version). The RPM and Debian package versions follow (`libopensps0`). The C
API version (`OPENSPS_API_VERSION_*`, `opensps_get_driver_version`, now 1.2) stays separate: it
describes the interface the bindings check, changes only when entry points do, and would
otherwise jump with every library release. The `FSPS Version: 3.2` line of the verbose banner
also stays: it names the FSPS version whose physics the code follows.

## Porting a patch from this repository after the rename

The trees match the libopensps root commit apart from `data/`, so a patch applies there and the
rename can be replayed:

```sh
cd libopensps
root=$(git rev-list --max-parents=0 HEAD)
git switch -c port-A "$root" && git am /path/to/0001-*.patch
tools/cutover/rename_to_opensps.sh . && git add -A && git commit -m "rename (port)"
git switch -c rename-only "$root" && tools/cutover/rename_to_opensps.sh . && git add -A && git commit -m "rename"
git switch main && git diff rename-only port-A | git apply -3
```
