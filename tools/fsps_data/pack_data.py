#!/usr/bin/env python3
"""Pack a standard-format data root into versioned archives for distribution.

Usage: pack_data.py ROOT OUTDIR --release RELEASE [--base-url URL] [--optional KIND/NAME ...]

Writes into OUTDIR one ``.tar.gz`` per library (``libraries/<kind>/<name>``) plus a base archive
(``fsps-data.json`` and the small tables in ``data/``), and ``fsps-data-<RELEASE>.json``, the
release manifest that fetch_data.py reads: every archive with its SHA-256, size and contents.
Archives hold paths relative to the data root, so unpacking all of them into one directory
gives the root back. Libraries listed with --optional (default: the C3K spectra, 0.6 GB) are
fetched only on request.

Every library is verified (sizes and checksums of its files) before it is packed.
"""
import argparse
import hashlib
import json
import os
import sys
import tarfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fsps_data import Library, read_root_index  # noqa: E402

DEFAULT_OPTIONAL = ["stellar_spectra/c3k_afe+0.0"]
# TODO(hosting): the archives are not published yet; set the URL of the release (e.g. a Zenodo
# record or a GitHub release) when the data are uploaded.
PLACEHOLDER_URL = "https://example.invalid/fsps-data/{release}/"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def add_tree(tar, root, rel):
    """Adds root/rel recursively with deterministic order and metadata."""
    full = os.path.join(root, rel)
    if os.path.isdir(full):
        for name in sorted(os.listdir(full)):
            add_tree(tar, root, os.path.join(rel, name))
        return
    info = tar.gettarinfo(full, arcname=rel)
    info.uid = info.gid = 0
    info.uname = info.gname = ""
    info.mtime = 0
    info.mode = 0o644
    with open(full, "rb") as f:
        tar.addfile(info, f)


def pack(root, outdir, fname, rels):
    path = os.path.join(outdir, fname)
    with tarfile.open(path, "w:gz", compresslevel=6) as tar:
        for rel in rels:
            add_tree(tar, root, rel)
    return {"file": fname, "bytes": os.path.getsize(path), "sha256": sha256(path), "contents": rels}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("outdir")
    ap.add_argument("--release", required=True)
    ap.add_argument("--base-url", default=None, help="where the archives will be hosted (placeholder if omitted)")
    ap.add_argument("--optional", nargs="*", default=DEFAULT_OPTIONAL)
    a = ap.parse_args(argv)

    idx = read_root_index(a.root)
    os.makedirs(a.outdir, exist_ok=True)
    archives = []
    base = ["fsps-data.json"] + [os.path.join("data", f) for f in sorted(os.listdir(os.path.join(a.root, "data")))]
    e = pack(a.root, a.outdir, f"fsps-data-{a.release}-base.tar.gz", base)
    e.update({"description": "root index and small text tables", "optional": False})
    archives.append(e)
    for lib in idx["libraries"]:
        key = f"{lib['kind']}/{lib['name']}"
        problems = Library.open(a.root, lib["kind"], lib["name"]).verify()
        if problems:
            sys.exit("not packing an inconsistent library:\n  " + "\n  ".join(problems))
        e = pack(a.root, a.outdir, f"fsps-data-{a.release}-{lib['kind']}-{lib['name']}.tar.gz", [lib["path"]])
        e.update({"library": key, "optional": key in a.optional})
        archives.append(e)
        print(f"  {e['file']}: {e['bytes'] / 1e6:.1f} MB", flush=True)
    rel = {
        "format": "fsps-data-release",
        "format_version": 1,
        "release": a.release,
        "data_format_version": idx["format_version"],
        "base_url": a.base_url or PLACEHOLDER_URL.format(release=a.release),
        "archives": archives,
    }
    out = os.path.join(a.outdir, f"fsps-data-{a.release}.json")
    with open(out, "w") as f:
        json.dump(rel, f, indent=1)
        f.write("\n")
    total = sum(x["bytes"] for x in archives)
    print(f"{len(archives)} archives, {total / 1e6:.0f} MB; release manifest {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
