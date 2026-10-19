#!/usr/bin/env python3
"""Download, verify and unpack an FSPS data release (doc/DATA_FORMAT.md, Section 6).

Usage:
    tools/fetch_data.py DEST [--release-manifest PATH_OR_URL] [--base-url URL]
                             [--with KIND/NAME ...] [--all] [--cache DIR] [--list]

Reads the release manifest (``fsps-data-<release>.json``, written by
tools/fsps_data/pack_data.py), downloads every archive that is not optional (plus those named
with --with, or all with --all) into a cache, checks each archive's SHA-256, unpacks it into
DEST, and finally checks the size and SHA-256 of every file the unpacked manifests list. Point
FSPS_DATA_HOME at DEST afterwards. Archives already in the cache with the right checksum are
not downloaded again (CI caches that directory).

Standard library only (no NumPy), Python >= 3.8.
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
import tarfile
import urllib.request

# TODO(hosting): the release manifest of the current data release; to be set when the archives
# are published (e.g. on Zenodo or as GitHub release assets).
DEFAULT_MANIFEST = "https://example.invalid/fsps-data/latest/fsps-data-latest.json"
CHUNK = 1 << 22


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(CHUNK), b""):
            h.update(b)
    return h.hexdigest()


def fetch(src, dest):
    """Copies a local path or downloads a URL to dest (through a temporary file)."""
    tmp = dest + ".part"
    if "://" not in src or src.startswith("file://"):
        shutil.copyfile(src.replace("file://", ""), tmp)
    else:
        with urllib.request.urlopen(src) as r, open(tmp, "wb") as f:
            shutil.copyfileobj(r, f, CHUNK)
    os.replace(tmp, dest)


def safe_extract(tar, dest):
    """Unpacks only regular files and directories below dest (no links, no absolute paths)."""
    root = os.path.realpath(dest)
    for m in tar.getmembers():
        target = os.path.realpath(os.path.join(dest, m.name))
        if not (m.isfile() or m.isdir()) or os.path.isabs(m.name) or \
                os.path.commonpath([root, target]) != root:
            raise ValueError(f"refusing archive member {m.name!r}")
    for m in tar.getmembers():
        tar.extract(m, dest)


def verify_root(dest):
    """Problems found in the files of every library listed in DEST/fsps-data.json."""
    problems = []
    with open(os.path.join(dest, "fsps-data.json")) as f:
        idx = json.load(f)
    for lib in idx["libraries"]:
        d = os.path.join(dest, lib["path"])
        mpath = os.path.join(d, "manifest.json")
        if not os.path.exists(mpath):
            continue  # an optional library that was not fetched
        with open(mpath) as f:
            man = json.load(f)
        for desc in list(man["axes"].values()) + list(man["arrays"].values()):
            if "file" not in desc:
                continue
            n = 1
            for s in desc["shape"]:
                n *= s
            nbytes = n * int(desc["dtype"][2:])
            p = os.path.join(d, desc["file"])
            if not os.path.exists(p) or os.path.getsize(p) != nbytes:
                problems.append(f"{p}: missing or wrong size")
            elif sha256(p) != desc["sha256"]:
                problems.append(f"{p}: checksum mismatch")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dest")
    ap.add_argument("--release-manifest", default=DEFAULT_MANIFEST)
    ap.add_argument("--base-url", default=None, help="override the archive location given in the manifest")
    ap.add_argument("--with", dest="extra", nargs="*", default=[], help="optional libraries to fetch too")
    ap.add_argument("--all", action="store_true", help="fetch every archive, optional ones included")
    ap.add_argument("--cache", default=os.path.join(os.path.expanduser("~"), ".cache", "fsps-data"))
    ap.add_argument("--list", action="store_true", help="list the archives of the release and exit")
    a = ap.parse_args(argv)

    if "example.invalid" in a.release_manifest:
        sys.exit("No data release is published yet (TODO: hosting). Pass --release-manifest PATH_OR_URL, "
                 "e.g. the fsps-data-<release>.json written by tools/fsps_data/pack_data.py.")
    os.makedirs(a.cache, exist_ok=True)
    mpath = os.path.join(a.cache, os.path.basename(a.release_manifest))
    fetch(a.release_manifest, mpath)
    with open(mpath) as f:
        rel = json.load(f)
    if rel.get("format") != "fsps-data-release" or rel.get("data_format_version") != 1:
        sys.exit(f"{a.release_manifest}: not a release of data format version 1")
    base = a.base_url or rel["base_url"]
    if "example.invalid" in base and a.base_url is None and "://" not in a.release_manifest:
        base = os.path.dirname(os.path.abspath(a.release_manifest))  # archives next to a local manifest
    if a.list:
        for e in rel["archives"]:
            print(f"{e['file']:60s} {e['bytes'] / 1e6:8.1f} MB{'  (optional)' if e['optional'] else ''}")
        return 0

    want = [e for e in rel["archives"] if a.all or not e["optional"] or e.get("library") in a.extra]
    os.makedirs(a.dest, exist_ok=True)
    for e in want:
        cached = os.path.join(a.cache, e["file"])
        if not (os.path.exists(cached) and os.path.getsize(cached) == e["bytes"] and sha256(cached) == e["sha256"]):
            src = base.rstrip("/") + "/" + e["file"] if "://" in base else os.path.join(base, e["file"])
            print(f"fetching {e['file']} ({e['bytes'] / 1e6:.1f} MB)", flush=True)
            fetch(src, cached)
            if sha256(cached) != e["sha256"]:
                os.remove(cached)
                sys.exit(f"{e['file']}: checksum mismatch after download")
        with tarfile.open(cached, "r:gz") as tar:
            safe_extract(tar, a.dest)
    problems = verify_root(a.dest)
    if problems:
        sys.exit("verification failed:\n  " + "\n  ".join(problems))
    print(f"release {rel['release']}: {len(want)} archives unpacked and verified in {a.dest}")
    print(f"export FSPS_DATA_HOME={os.path.abspath(a.dest)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
