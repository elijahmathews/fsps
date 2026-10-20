# Cutover day checklist

The steps in order, with commands. `FORK` is a clone of elijahmathews/fsps, `NEW` the clone of
libopensps where the import is built. Scripts and decisions: [README.md](README.md). Every `gh`
command names its repository with `-R` / `--repo`.

```sh
FORK=~/Code/fsps-fork-clone
NEW=~/Code/libopensps
export OPENSPS_DATA_HOME=~/opensps-data
```

Decided: version 0.1.0 (`libopensps.so.0`, the rename script's default; the C API version stays
1.2); the regression references are left out of the import (the default) for the data release;
the two Google Code era handles stay handles. Open: whether the internal Fortran names follow
(README.md, "Renaming the internal names").

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

## 3. Build the import commit locally

```sh
git clone https://github.com/elijahmathews/libopensps.git $NEW
$FORK/tools/cutover/make_import_commit.sh --repo $FORK $CUT $NEW    # branch `import`, one commit
cd $NEW && git log --format='%h %an %s' -5 && git show --stat HEAD | head -30 && cat AUTHORS
```

The commit is `feat: import elijahmathews/fsps at commit $CUT` on top of `main`, unsigned,
authored by your git identity, without trailers. `main`'s own README.md, LICENSE.md, AGENTS.md
(and CLAUDE.md) are kept; `proposals/README.md` has README sections to merge by hand later. The
regression references left out are in `$NEW.references`; check them against data-v2 and run the
suites from the fetched release:

```sh
python3 tools/fetch_data.py $OPENSPS_DATA_HOME --release data-v2 --all
(cd $NEW.references && for f in $(find . -type f ! -name SHA256SUMS); do
   cmp -s "$f" "$(find $OPENSPS_DATA_HOME/libraries/references -name "$(basename "$f")" | head -1)" || echo "DIFFERS: $f"; done)
meson setup build --buildtype=release && meson setup build32 --buildtype=release -Dkernel_precision=single
meson compile -C build && meson compile -C build32
meson test -C build --suite unit --suite regression && meson test -C build32 --suite unit --suite regression
```

## 4. Land the import on main

First, so that CI on the pull request fetches the data and the references:

```sh
gh variable set FSPS_DATA_RELEASE -R elijahmathews/libopensps --body data-v2
```

**Recommended: a pull request, squash-merged** (one commit stays one commit, with its message):

```sh
cd $NEW && git push -u origin import
gh pr create -R elijahmathews/libopensps --base main --head import   --title "feat: import elijahmathews/fsps at commit $CUT" --body "$(git log -1 --format=%b)"
gh pr checks -R elijahmathews/libopensps import --watch
gh pr merge -R elijahmathews/libopensps import --squash   --subject "feat: import elijahmathews/fsps at commit $CUT" --body "$(git log -1 --format=%b)"
```

(or merge in the web interface: check that the squash message is the commit's.)

**Alternative: a direct push** by you (bypasses the protection and CI before landing; the
commit lands exactly as built):

```sh
cd $NEW && git push origin import:main
```

## 5. The rename (its own pull request)

```sh
cd $NEW && git fetch && git switch -c refactor/rename-opensps origin/main
tools/cutover/rename_to_opensps.sh .                  # version 0.1.0, libopensps.so.0
tools/cutover/rename_to_opensps.sh --check-only . && git diff --stat
meson setup build --wipe && meson compile -C build && meson test -C build --suite unit --suite regression
readelf -d build/libopensps.so | grep SONAME          # libopensps.so.0
nm -D --defined-only build/libopensps.so | awk '{print $3}' | grep -E '^fsps\w*[^_]$'   # no C functions: expect nothing
git add -A && git -c commit.gpgsign=false commit -F - <<'MSG'
refactor!: rename the public C API and library to opensps

The library is libopensps (version 0.1.0, soname libopensps.so.0), the header
opensps.h, the pkg-config name opensps; C functions and macros carry the
prefixes opensps_ and OPENSPS_. Internal Fortran names are unchanged.

BREAKING CHANGE: libfsps, fsps.h, pkg-config fsps and the fsps_* / FSPS_* C
names are now libopensps, opensps.h, opensps and opensps_* / OPENSPS_*. A
deprecated fsps.h maps the old names; link against libopensps.
MSG
git push -u origin refactor/rename-opensps
gh pr create -R elijahmathews/libopensps --base main --head refactor/rename-opensps --fill
gh pr checks -R elijahmathews/libopensps refactor/rename-opensps --watch
```

Squash-merge once CI passes (the pull request title is the commit subject: keep it the
Conventional one above, with the `BREAKING CHANGE:` footer in the body).

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
