# Cutover day checklist

The steps in order, with commands. `FORK` is a clone of elijahmathews/fsps, `NEW` the directory
where libopensps is assembled. Scripts and decisions: [README.md](README.md). Every `gh`
command names its repository with `-R` / `--repo`.

```sh
FORK=~/Code/fsps-fork-clone
NEW=~/Code/libopensps
export OPENSPS_DATA_HOME=~/opensps-data
```

Decided: version 0.1.0 (`libopensps.so.0`, the rename script's default; the C API version stays
1.2); the regression references leave the snapshot (the default) for the data release; the two
Google Code era handles stay handles.

## 0. Before the day

- [ ] **Regression references in the data release** (branch `data/references`): the tests
      resolve `sps_ref_*.bin` and `param_ref/` through the data root, falling back to
      `tests/data/`. It must be merged before the cutover commit, since the snapshot leaves
      `tests/data/sps_ref_*.bin` and `tests/data/param_ref/` out.
- [ ] **Publish data-v2** (the release that contains the references) to
      [elijahmathews/opensps-data](https://github.com/elijahmathews/opensps-data). This pushes to
      a public repository: only with your explicit go-ahead.

  ```sh
  cat ~/Code/fsps-data-release/data-v2/UPLOAD.txt                     # what will be uploaded
  $FORK/tools/fsps_data/upload_release.sh ~/Code/fsps-data-release/data-v2 data-v2   # asks before uploading
  python3 $FORK/tools/fetch_data.py /tmp/check-v2 --release data-v2 --all --cache /tmp/check-v2-cache
  ```

- [ ] **Default release** in `tools/fetch_data.py` (`DEFAULT_RELEASE`) set to `data-v2`, so a
      plain `fetch_data.py ~/opensps-data` includes the references.
- [ ] **CI of this repository on data-v2:** the repository variable `FSPS_DATA_RELEASE` (read by
      `ci.yml`) set to `data-v2`, and one CI run on develop that fetches the release.

  ```sh
  gh variable set FSPS_DATA_RELEASE -R elijahmathews/fsps --body data-v2
  gh workflow run ci.yml -R elijahmathews/fsps --ref develop
  ```

## 1. Merge into develop (order decided by you)

Whatever is ready, including `data/references` and this pull request (`tools/cutover` must be in
the cutover commit, so that the rename script travels with the tree). Work not merged moves later
as patches (README.md, "Porting a patch").

## 2. The cutover commit

```sh
cd $FORK && git switch develop && git pull --ff-only
gh run list -R elijahmathews/fsps --branch develop --limit 5         # CI green on this commit
gh workflow run nightly.yml -R elijahmathews/fsps --ref develop      # and the nightly suites
git tag -a fsps-final -m "Last commit before the move to elijahmathews/libopensps"
git push origin fsps-final
CUT=$(git rev-parse fsps-final^{commit})
```

## 3. Build libopensps locally

```sh
$FORK/tools/cutover/make_initial_tree.sh --repo $FORK $CUT $NEW
cd $NEW && git log --stat -1 | head -20 && cat AUTHORS && head -6 LICENSE.md
du -sh .git                                   # a few MB; the references were copied to $NEW.references
```

The commit is "libopensps: the state of elijahmathews/fsps at commit $CUT", unsigned, authored
by your git identity, without trailers (`--trailer` adds some if wanted). Check that the
references left out are those of data-v2, then build and test from the fetched release:

```sh
python3 tools/fetch_data.py $OPENSPS_DATA_HOME --release data-v2 --all
# the references library of the release against $NEW.references/SHA256SUMS (same files, same sums)
meson setup build --buildtype=release && meson setup build32 --buildtype=release -Dkernel_precision=single
meson compile -C build && meson compile -C build32
meson test -C build --suite unit --suite regression && meson test -C build32 --suite unit --suite regression
```

## 4. Push to libopensps

```sh
cd $NEW
git remote add origin git@github.com:elijahmathews/libopensps.git
git push -u origin main
gh variable set FSPS_DATA_RELEASE -R elijahmathews/libopensps --body data-v2   # CI fetches the references with the data
```

## 5. The rename (first pull request in libopensps)

```sh
cd $NEW && git switch -c rename/opensps
tools/cutover/rename_to_opensps.sh .                  # version 0.1.0, libopensps.so.0
tools/cutover/rename_to_opensps.sh --check-only . && git diff --stat
meson setup build --wipe && meson compile -C build && meson test -C build --suite unit --suite regression
readelf -d build/libopensps.so | grep SONAME          # libopensps.so.0
nm -D --defined-only build/libopensps.so | awk '{print $3}' | grep -i '^fsps'   # expect nothing
git add -A && git -c commit.gpgsign=false commit -m "Rename the public surface to opensps; version 0.1.0"
git push -u origin rename/opensps
gh pr create -R elijahmathews/libopensps --base main --head rename/opensps --fill
gh pr checks -R elijahmathews/libopensps rename/opensps --watch
```

Review the installed `README.md`, `CITATION.cff` (ORCID, version, DOI are TODO) and
`CONTRIBUTING.md`; merge once CI passes.

## 6. Repository settings (github.com/elijahmathews/libopensps/settings)

- [ ] Description ("Stellar population synthesis library derived from FSPS"), topics
      (`stellar-population-synthesis`, `astronomy`, `fortran`, `sed`).
- [ ] Issues on; wiki and projects off.
- [ ] Default branch `main`; branch protection / ruleset on `main`: pull requests required, the
      CI checks required, no force pushes.
- [ ] Actions allowed; variable `FSPS_DATA_RELEASE` = `data-v2` (step 4); secret
      `CODECOV_TOKEN`; anything the nightly or data jobs use.
- [ ] Optionally connect the repository to an archive service that mints a DOI per release, and
      fill the DOI in `CITATION.cff`.

```sh
gh repo edit elijahmathews/libopensps --enable-issues --enable-wiki=false --enable-projects=false \
  --description "Stellar population synthesis library derived from FSPS"
gh secret set CODECOV_TOKEN -R elijahmathews/libopensps
```

## 7. First release, then the Python bindings

**Tag v0.1.0 only after the rename pull request has merged and CI on `main` is green** for the
merge commit:

```sh
gh run list -R elijahmathews/libopensps --branch main --limit 3      # green on the merge commit
cd $NEW && git switch main && git pull --ff-only
git tag -a v0.1.0 -m "libopensps 0.1.0" && git push origin v0.1.0
gh release create v0.1.0 -R elijahmathews/libopensps --generate-notes
```

In opensps-python (one pull request), pinned to that tag:

- `src/opensps/_capi.pxd`: `fsps.h` -> `opensps.h`, every `"fsps_..."` -> `"opensps_..."`,
  `"FSPS_..."` -> `"OPENSPS_..."` (the only file naming C symbols);
- `subprojects/fsps.wrap` -> `subprojects/opensps.wrap`: `url =
  https://github.com/elijahmathews/libopensps.git`, `revision = v0.1.0`, `[provide] opensps =
  opensps_dep`;
- `meson.build`: `fallback : ['opensps', 'opensps_dep']`; docs: `-Dfsps:` -> `-Dopensps:`,
  `subprojects/fsps` -> `subprojects/opensps`, and the note about the subproject's data download
  (it no longer carries data).

## 8. Retire the fork

```sh
cd $FORK && git switch -c docs/moved develop
# README.md: a note at the top: "Development continues at https://github.com/elijahmathews/libopensps"
git -c commit.gpgsign=false commit -am "README: development moved to libopensps" && git push -u origin docs/moved
gh pr create -R elijahmathews/fsps --base develop --head docs/moved --fill   # merge it, then:
gh repo archive elijahmathews/fsps --yes
```

The tag `fsps-final` keeps the full history and the data reachable.
