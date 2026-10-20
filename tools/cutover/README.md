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
`.mod` files, `fsps-data/`) can get in. Of the tracked files it leaves out:

- **`data/`** (1.4 GB): the libraries, distributed as versioned archives (PR #21,
  `tools/fetch_data.py`).

Files above 1 MB that remain (the script lists them on every run):

| File | Size | Decision |
|---|---|---|
| `tests/data/sps_ref_*.bin` (5 files) | 52 MB | **kept** by default; `--exclude-refs DIR` moves them out |
| `tests/data/fsps-test-data/.../miles/flux.bin` (after PR #21) | 2.4 MB | kept: the test library the unit suite runs on |

Why the regression references stay in the first commit: they are the test oracles of the
regression suite (meson registers them by path), so removing them breaks `meson test --suite
regression` until a follow-up teaches meson to find them elsewhere; and the roadmap ports
unfinished work with `git format-patch` / `git am`, which relies on the trees being identical
apart from `data/`. With them, the repository is 47 MB of git objects (49 MB with PR #21),
within the roadmap's ~50 MB goal. The cost comes later: each regeneration of the references
adds about 50 MB to libopensps's history for good. Recommendation: move them out (to the data
release, or a separate reference archive fetched like the data) in libopensps **before** they
are first regenerated. Without them the repository would be 7 MB (10 MB with PR #21).
`--exclude-refs DIR` does it at cutover instead and writes the files with `SHA256SUMS` to `DIR`.

Everything else is small (< 1 MB): `doc/MANUAL.pdf` (0.16 MB, the FSPS manual), the
parameterized references `tests/data/param_ref/*.opsf` (7.1 MB together, each below 1 MB), editor and agent
files (`.vscode/`, `.fortls`, `CLAUDE.md`, `GEMINI.md`), all kept so that the snapshot is the
fork's tree.

## AUTHORS and the license

The history that records authorship is not carried over, so the snapshot adds `AUTHORS`,
generated from the history of `<commit>` (`git log --use-mailmap`, ordered by number of
commits): FSPS contributors are those with commits up to the fork point (`--fsps-base`, the
merge base with `upstream/master`, `82a873508d50`), fork contributors those after it. Identities
are merged only where they are certainly one person (`authors.mailmap`: Charlie Conroy under six
addresses, Benjamin Johnson, Nell Byler, Elijah Mathews under two each). Two Google Code era
committers appear only as addresses (`tomorrowsville`, 2011-2013; `avillaum`, 2014) and are
listed as such; if they can be identified, add them to `authors.mailmap` and rerun. Bots and AI
agents are not listed as authors; a sentence notes that the fork was written with AI agents.

`LICENSE.md` keeps the FSPS notice unchanged and gains one line after it:
`Copyright © 2026 Elijah Mathews & contributors (the modernized fork, libopensps).`

## Rename scope

**Recommendation: rename the public surface, keep the internal Fortran names.** Implemented:

| Renamed | From | To |
|---|---|---|
| Library and soname | `libfsps.so.3` | `libopensps.so.3` (`--version X.Y.Z` sets version and soversion) |
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
  directory names (`share/fsps`, `~/.local/share/fsps`): existing setups, CI and the Python
  bindings' data search keep working. An `OPENSPS_DATA_HOME` alias would be a small code
  change, better made deliberately than by a rename script.
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
