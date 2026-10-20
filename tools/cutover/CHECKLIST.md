# Cutover day checklist

The steps in order, with commands. `FORK` is a clone of elijahmathews/fsps, `NEW` the directory
where libopensps is assembled. Scripts and decisions: [README.md](README.md).

```sh
FORK=~/Code/fsps-fork-clone
NEW=~/Code/libopensps
```

## 0. Before the day (blockers found in the dry run)

- [ ] **Unit suite without data.** `unit_io` ("Nebular Grid Interpolation" in
      `tests/io/test_fsps_io.f90`) reads `data/nebular/ZAU_WD_mist.cont` from `FSPS_HOME`,
      `MESON_SOURCE_ROOT` or the working directory, so the unit suite fails in a tree without
      `data/` (also with PR #21 at 8cc4455: 22 of 23 unit tests pass). Point it at the test
      library or skip it when the file is absent.
- [ ] **Regression tests and the data location.** `meson.build` sets `FSPS_DATA_HOME` to the
      source root for the regression suite, and `test_physics_oracles` reads
      `$FSPS_DATA_HOME/data/emlines_info.dat` directly. Without `data/` in the tree the
      regression suite only works because `SPS_HOME` overrides the library's search, and
      `physics_oracles` still fails ("H-beta found in the line list"). Let the environment (or a
      `-Ddata_home=` option) set the data root and have the test ask the library for its path.
      The C3K regression is registered only if `data/spectra/C3K` exists in the source tree.
- [ ] **CI without `data/`.** `ci.yml` checks out `data/` (sparse checkout without C3K); in
      libopensps the jobs must fetch the release (`doc/DATA_FORMAT.md`, section 6). Same for the
      packaging workflow's `%check` / `dh_auto_test`, which run the regression suite.
- [ ] **Data release published** (`data-v1`) and the TODO URLs in `tools/fetch_data.py` /
      `pack_data.py` set; the README template's fetch command updated.
- [ ] Optional: identify the Google Code committers `tomorrowsville` and `avillaum` (add them
      to `authors.mailmap`); decide the first libopensps version (`--version`, e.g. `0.1.0`,
      soname `libopensps.so.0`; default keeps 3.2.0 / `.so.3`); decide on `--exclude-refs`
      (README.md, "What the first commit leaves out").

## 1. Merge into develop (order decided by you)

Open, ready pull requests (as of 2026-10-20): #20 C API for bindings, #24 batch young/old (on
#20), #23 tabular SFH mass formed, #21 standard data format, #19 roadmap, #22 scenario as batch
(opt-in), #25 CI retry, #13 / #16 OpenMP mapping for LLVM flang, and this one (tools/cutover).
This pull request must be in the cutover commit, so that the rename script travels with the
tree. Anything not merged moves later as patches (README.md, "Porting a patch").

## 2. The cutover commit

```sh
cd $FORK && git switch develop && git pull --ff-only
# wait for a green CI run on this commit, including the nightly suites (workflow_dispatch)
gh run list --branch develop --limit 5
git tag -a fsps-final -m "Last commit before the move to elijahmathews/libopensps"
git push origin fsps-final
CUT=$(git rev-parse fsps-final^{commit})
```

## 3. Build libopensps locally

```sh
$FORK/tools/cutover/make_initial_tree.sh --repo $FORK $CUT $NEW   # add --exclude-refs DIR if decided
cd $NEW && git log --stat -1 | head -20 && cat AUTHORS && head -6 LICENSE.md
du -sh .git                                                        # ~47-49 MB (7-10 MB without refs)
```

The commit is "libopensps: the state of elijahmathews/fsps at commit $CUT", unsigned, authored
by your git identity, without trailers (`--trailer` adds some if wanted).

Check it: build and test from fetched data.

```sh
python3 tools/fetch_data.py ~/opensps-data && export FSPS_DATA_HOME=~/opensps-data
meson setup build --buildtype=release && meson setup build32 --buildtype=release -Dkernel_precision=single
meson compile -C build && meson compile -C build32
meson test -C build --suite unit --suite regression && meson test -C build32 --suite unit --suite regression
```

## 4. Push to libopensps

```sh
cd $NEW
git remote add origin git@github.com:elijahmathews/libopensps.git
git push -u origin main
```

## 5. The rename (first pull request in libopensps)

```sh
cd $NEW && git switch -c rename/opensps
tools/cutover/rename_to_opensps.sh .            # add --version X.Y.Z if decided
git diff --stat && tools/cutover/rename_to_opensps.sh --check-only .
meson setup build --wipe && meson compile -C build && meson test -C build --suite unit --suite regression
nm -D --defined-only build/libopensps.so | awk '{print $3}' | grep -i '^fsps'   # expect nothing
git add -A && git -c commit.gpgsign=false commit -m "Rename the public surface to opensps"
git push -u origin rename/opensps && gh pr create --base main --fill
```

Review the installed `README.md`, `CITATION.cff` (ORCID, version, DOI are TODO) and
`CONTRIBUTING.md`; merge after CI passes.

## 6. Repository settings (github.com/elijahmathews/libopensps/settings)

- [ ] Description ("Stellar population synthesis library derived from FSPS"), website, topics
      (`stellar-population-synthesis`, `astronomy`, `fortran`, `sed`).
- [ ] Issues on; discussions as wanted; wiki and projects off.
- [ ] Default branch `main`; branch protection / ruleset on `main`: pull requests required,
      the CI checks required, no force pushes, linear history if wanted.
- [ ] Actions: allow; secrets `CODECOV_TOKEN` (codecov.yml) and anything the nightly GPU jobs
      or data upload use; the coverage badge points at libopensps.
- [ ] Optionally connect the repository to an archive service that mints a DOI per release, and
      fill the DOI in `CITATION.cff`.

```sh
gh repo edit elijahmathews/libopensps --enable-issues --enable-wiki=false \
  --description "Stellar population synthesis library derived from FSPS"
gh secret set CODECOV_TOKEN -R elijahmathews/libopensps
```

## 7. First release and the Python bindings

```sh
cd $NEW && git switch main && git pull
git tag -a v0.1.0 -m "libopensps 0.1.0" && git push origin v0.1.0   # the version chosen above
gh release create v0.1.0 -R elijahmathews/libopensps --generate-notes
```

In opensps-python (one pull request):

- `src/opensps/_capi.pxd`: `fsps.h` -> `opensps.h`, every `"fsps_..."` -> `"opensps_..."`,
  `"FSPS_..."` -> `"OPENSPS_..."` (the only file naming C symbols);
- `subprojects/fsps.wrap` -> `subprojects/opensps.wrap`: `url =
  https://github.com/elijahmathews/libopensps.git`, `revision = v0.1.0`, `[provide] opensps =
  opensps_dep`;
- `meson.build`: `fallback : ['opensps', 'opensps_dep']`; docs: `-Dfsps:` -> `-Dopensps:`,
  `subprojects/fsps` -> `subprojects/opensps`, the size note about the 1.3 GB data fetch;
- `_config.py` needs no change (the environment variables are unchanged).

## 8. Retire the fork

```sh
cd $FORK && git switch -c docs/moved develop
# README.md: a note at the top: "Development continues at https://github.com/elijahmathews/libopensps"
git -c commit.gpgsign=false commit -am "README: development moved to libopensps" && git push -u origin docs/moved
gh pr create --base develop --fill    # merge it, then:
gh repo archive elijahmathews/fsps
```

The tag `fsps-final` keeps the full history and the data reachable.
