#!/usr/bin/env python3
"""Convert legacy FSPS spectral binaries into an FSDS-style HDF5 file.

Usage:
  python tools/convert_spectra_to_fsds.py \
      --sps-home /path/to/fsps \
      --spec-lib miles \
    --isoc-lib mist \
      --out /path/to/output.h5
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


QPAH_ARR_DL07 = np.array([0.47, 1.12, 1.77, 2.50, 3.19, 3.90, 4.58], dtype=np.float64)
UMIN_ARR_DL07 = np.array(
    [
        0.1,
        0.15,
        0.2,
        0.3,
        0.4,
        0.5,
        0.7,
        0.8,
        1.0,
        1.2,
        1.5,
        2.0,
        2.5,
        3.0,
        4.0,
        5.0,
        7.0,
        8.0,
        12.0,
        15.0,
        20.0,
        25.0,
    ],
    dtype=np.float64,
)
QPAH_ARR_THEMIS = (100.0 / 2.2) * np.array(
    [0.02, 0.06, 0.10, 0.14, 0.17, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40],
    dtype=np.float64,
)
UMIN_ARR_THEMIS = np.array(
    [
        0.1,
        0.12,
        0.15,
        0.17,
        0.2,
        0.25,
        0.3,
        0.35,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        1.0,
        1.2,
        1.5,
        1.7,
        2.0,
        2.5,
        3.0,
        3.5,
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
        10.0,
        12.0,
        15.0,
        17.0,
        20.0,
        25.0,
        30.0,
        35.0,
        40.0,
        50.0,
        80.0,
    ],
    dtype=np.float64,
)


def _spectra_dir(sps_home: Path, spec_lib: str) -> Path:
    spec = spec_lib.lower()
    if spec == "bpass":
        return sps_home / "data" / "isochrones" / "BPASS"
    if spec == "miles":
        return sps_home / "data" / "spectra" / "MILES"
    if spec == "basel":
        return sps_home / "data" / "spectra" / "BaSeL3.1"
    if spec == "ckc14":
        return sps_home / "data" / "spectra" / "CKC14"
    if spec.startswith("c3k"):
        return sps_home / "data" / "spectra" / "C3K"
    raise ValueError(f"Unsupported spectral library: {spec_lib}")


def _lambda_filename(spec_lib: str) -> str:
    spec = spec_lib.lower()
    if spec == "bpass":
        return "bpass.lambda"
    if spec == "miles":
        return "miles.lambda"
    if spec == "basel":
        return "basel.lambda"
    if spec == "ckc14":
        return "ckc14.lambda"
    if spec.startswith("c3k"):
        return f"{spec_lib}.lambda"
    raise ValueError(f"Unsupported spectral library: {spec_lib}")


def _binary_filename(spec_lib: str, z: float, spectra_dir: Path) -> Path:
    zstr = f"{z:0.4f}"
    spec = spec_lib.lower()

    if spec == "bpass":
        return spectra_dir / "bpass_v2.2_salpeter100.ssp.bin"

    if spec == "miles":
        return spectra_dir / f"imiles_z{zstr}.spectra.bin"

    if spec == "basel":
        # BASEL_STR is build/runtime dependent in FSPS. Try direct pattern first,
        # then fall back to a wildcard match.
        direct = spectra_dir / f"basel_z{zstr}.spectra.bin"
        if direct.exists():
            return direct
        matches = sorted(spectra_dir.glob(f"basel_*_z{zstr}.spectra.bin"))
        if matches:
            return matches[0]
        return direct

    if spec == "ckc14" or spec.startswith("c3k"):
        return spectra_dir / f"{spec_lib}_z{zstr}.spectra.bin"

    raise ValueError(f"Unsupported spectral library: {spec_lib}")


def _isoc_dir(sps_home: Path, isoc_lib: str) -> Path:
    lib = isoc_lib.lower()
    if lib == "bpss":
        return sps_home / "data" / "isochrones" / "BPASS"
    if lib == "mist":
        return sps_home / "data" / "isochrones" / "MIST"
    if lib == "pdva":
        return sps_home / "data" / "isochrones" / "Padova" / "Padova2007"
    if lib == "prsc":
        return sps_home / "data" / "isochrones" / "PARSEC"
    if lib == "bsti":
        return sps_home / "data" / "isochrones" / "BaSTI"
    if lib == "gnva":
        return sps_home / "data" / "isochrones" / "Geneva"
    raise ValueError(f"Unsupported isochrone library: {isoc_lib}")


def _isoc_zsol(isoc_lib: str) -> float:
    lib = isoc_lib.lower()
    if lib == "bpss":
        return 0.020
    if lib == "mist":
        return 0.0142
    if lib == "pdva":
        return 0.019
    if lib == "prsc":
        return 0.01524
    if lib == "bsti":
        return 0.020
    if lib == "gnva":
        return 0.020
    raise ValueError(f"Unsupported isochrone library: {isoc_lib}")


def _read_axis(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Missing axis file: {path}")
    return np.loadtxt(path, dtype=np.float64, ndmin=1)


def _read_isoc_zlegend(isoc_dir: Path, isoc_lib: str) -> np.ndarray:
    path = isoc_dir / "zlegend.dat"
    if not path.exists():
        raise FileNotFoundError(f"Missing isochrone zlegend file: {path}")

    lib = isoc_lib.lower()
    if lib != "mist":
        return np.loadtxt(path, dtype=np.float64, ndmin=1)

    zsol = _isoc_zsol(lib)
    vals: list[float] = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        tag = s[0].lower()
        val = float(s[1:])
        if tag == "m":
            vals.append((10.0 ** (-val)) * zsol)
        else:
            vals.append((10.0 ** (val)) * zsol)

    if not vals:
        raise ValueError(f"No entries found in {path}")
    return np.asarray(vals, dtype=np.float64)


def _iso_filename(isoc_dir: Path, isoc_lib: str, z_val: float) -> Path:
    lib = isoc_lib.lower()
    if lib == "mist":
        val_log = np.log10(z_val / 0.0142)
        z_str = f"m{abs(val_log):0.2f}" if val_log < -0.001 else f"p{abs(val_log):0.2f}"
    else:
        z_str = f"{z_val:0.4f}"
    return isoc_dir / f"isoc_z{z_str}.dat"


def _parse_iso_file(
    path: Path, is_mist: bool
) -> list[list[tuple[float, float, float, float, float, float, float, float, float]]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing isochrone file: {path}")

    tracks: list[
        list[tuple[float, float, float, float, float, float, float, float, float]]
    ] = []
    start_new_track = True

    for raw in path.read_text().splitlines():
        line = raw.lstrip()
        if not line:
            continue

        if line.startswith("#"):
            start_new_track = True
            continue

        parts = line.split()
        try:
            if is_mist:
                if len(parts) < 9:
                    start_new_track = True
                    continue
                logage, mini, mact, logl, logt, logg, ffco, phase, lmdot = map(
                    float, parts[:9]
                )
            else:
                if len(parts) < 8:
                    start_new_track = True
                    continue
                logage, mini, mact, logl, logt, logg, ffco, phase = map(
                    float, parts[:8]
                )
                lmdot = -99.0
        except ValueError:
            start_new_track = True
            continue

        if start_new_track:
            tracks.append([])
            start_new_track = False

        tracks[-1].append((logage, mini, mact, logl, logt, logg, ffco, phase, lmdot))

    return [t for t in tracks if t]


def _read_legacy_cube(
    bin_path: Path, n_lam: int, n_logt: int, n_logg: int
) -> np.ndarray:
    if not bin_path.exists():
        raise FileNotFoundError(f"Missing binary spectral file: {bin_path}")

    expected = n_lam * n_logt * n_logg
    data = np.fromfile(bin_path, dtype=np.float32)
    if data.size != expected:
        raise ValueError(
            f"Unexpected data size in {bin_path}: got {data.size}, expected {expected}"
        )

    # Required by task: Fortran direct-access record interpreted in Python as
    # (n_logg, n_logt, n_lam).
    return data.reshape((n_logg, n_logt, n_lam))


def _read_bpass_mass_table(mass_path: Path, n_z: int) -> tuple[np.ndarray, np.ndarray]:
    if not mass_path.exists():
        raise FileNotFoundError(f"Missing BPASS mass file: {mass_path}")

    raw = np.fromstring(mass_path.read_text(), sep=" ", dtype=np.float64)
    n_cols = n_z + 1
    if raw.size == 0 or raw.size % n_cols != 0:
        raise ValueError(
            f"Unexpected BPASS mass table size in {mass_path}: got {raw.size}, expected multiple of {n_cols}"
        )

    arr = raw.reshape((-1, n_cols))
    time_full = arr[:, 0].astype(np.float64, copy=False)
    mass_ssp = arr[:, 1:].astype(np.float64, copy=False)
    return time_full, mass_ssp


def _read_bpass_cube(bin_path: Path, n_lam: int, n_t: int, n_z: int) -> np.ndarray:
    if not bin_path.exists():
        raise FileNotFoundError(f"Missing BPASS spectral binary file: {bin_path}")

    expected = n_lam * n_t * n_z
    data = np.fromfile(bin_path, dtype=np.float64)
    if data.size != expected:
        raise ValueError(
            f"Unexpected BPASS binary size in {bin_path}: got {data.size}, expected {expected}"
        )

    return data.reshape((n_lam, n_t, n_z), order="F")


def _write_bpss_isochrones(
    h5: h5py.File, isoc_lib: str, axis_z: np.ndarray, time_full: np.ndarray, mass_ssp: np.ndarray
) -> None:
    n_t = int(time_full.size)
    n_z = int(axis_z.size)
    if mass_ssp.shape != (n_t, n_z):
        raise ValueError(
            f"BPASS mass table shape mismatch: got {mass_ssp.shape}, expected {(n_t, n_z)}"
        )

    nmass = np.ones((n_z, n_t), dtype=np.int32)
    timestep_logyr = np.repeat(time_full[np.newaxis, :], n_z, axis=0).astype(
        np.float64, copy=False
    )

    missing64 = np.float64(-1.0e30)
    mini = np.full((n_z, n_t, 1), missing64, dtype=np.float64)
    mact = np.full((n_z, n_t, 1), missing64, dtype=np.float64)
    logl = np.full((n_z, n_t, 1), 0.0, dtype=np.float64)
    logt = np.full((n_z, n_t, 1), 0.0, dtype=np.float64)
    logg = np.full((n_z, n_t, 1), 0.0, dtype=np.float64)
    phase = np.full((n_z, n_t, 1), 0.0, dtype=np.float64)
    ffco = np.full((n_z, n_t, 1), 0.0, dtype=np.float64)
    lmdot = np.full((n_z, n_t, 1), -99.0, dtype=np.float64)

    mini[:, :, 0] = mass_ssp.T
    mact[:, :, 0] = mass_ssp.T

    tracks_grp = h5.require_group(f"/libraries/isochrones/{isoc_lib}/tracks")

    d_nmass = tracks_grp.create_dataset("nmass", data=nmass, dtype=np.int32)
    d_nmass.attrs["role"] = "isoc_nmass"
    d_nmass.attrs["dims_csv"] = "nt,nz"
    d_nmass.attrs["representation"] = "dense_nd"

    d_tstep = tracks_grp.create_dataset(
        "timestep", data=timestep_logyr, dtype=np.float64
    )
    d_tstep.attrs["role"] = "isoc_timestep"
    d_tstep.attrs["dims_csv"] = "nt,nz"
    d_tstep.attrs["representation"] = "dense_nd"

    def _write_3d(name: str, arr: np.ndarray) -> None:
        ds = tracks_grp.create_dataset(name, data=arr, dtype=np.float64)
        ds.attrs["role"] = f"isoc_{name}"
        ds.attrs["dims_csv"] = "nm,nt,nz"
        ds.attrs["representation"] = "dense_nd"

    _write_3d("mini", mini)
    _write_3d("mact", mact)
    _write_3d("logl", logl)
    _write_3d("logt", logt)
    _write_3d("logg", logg)
    _write_3d("phase", phase)
    _write_3d("ffco", ffco)
    _write_3d("lmdot", lmdot)


def _write_isochrones(h5: h5py.File, sps_home: Path, isoc_lib: str) -> None:
    isoc_dir = _isoc_dir(sps_home, isoc_lib)
    zlegend = _read_isoc_zlegend(isoc_dir, isoc_lib)

    if isoc_lib.lower() == "bpss":
        time_full, mass_ssp = _read_bpass_mass_table(isoc_dir / "bpass.mass", int(zlegend.size))
        _write_bpss_isochrones(h5, isoc_lib, zlegend, time_full, mass_ssp)
        return

    n_z = int(zlegend.size)
    is_mist = isoc_lib.lower() == "mist"

    all_tracks: list[
        list[list[tuple[float, float, float, float, float, float, float, float, float]]]
    ] = []
    max_nt = 0
    max_nm = 0

    for iz in range(n_z):
        z_val = float(zlegend[iz])
        path = _iso_filename(isoc_dir, isoc_lib, z_val)
        tracks = _parse_iso_file(path, is_mist)
        all_tracks.append(tracks)
        max_nt = max(max_nt, len(tracks))
        if tracks:
            max_nm = max(max_nm, max(len(t) for t in tracks))

    if max_nt <= 0 or max_nm <= 0:
        raise ValueError(f"No isochrone tracks parsed for {isoc_lib}")

    missing64 = np.float64(-1.0e30)
    nmass = np.zeros((n_z, max_nt), dtype=np.int32)
    timestep_logyr = np.full((n_z, max_nt), -1.0e30, dtype=np.float64)

    mini = np.full((n_z, max_nt, max_nm), missing64, dtype=np.float64)
    mact = np.full((n_z, max_nt, max_nm), missing64, dtype=np.float64)
    logl = np.full((n_z, max_nt, max_nm), missing64, dtype=np.float64)
    logt = np.full((n_z, max_nt, max_nm), missing64, dtype=np.float64)
    logg = np.full((n_z, max_nt, max_nm), missing64, dtype=np.float64)
    phase = np.full((n_z, max_nt, max_nm), missing64, dtype=np.float64)
    ffco = np.full((n_z, max_nt, max_nm), missing64, dtype=np.float64)
    lmdot = np.full((n_z, max_nt, max_nm), -99.0, dtype=np.float64)

    for iz, tracks in enumerate(all_tracks):
        for it, track in enumerate(tracks):
            if not track:
                continue
            nm = len(track)
            nmass[iz, it] = nm
            timestep_logyr[iz, it] = track[0][0]
            for im, row in enumerate(track):
                _, v_mini, v_mact, v_logl, v_logt, v_logg, v_ffco, v_phase, v_lmdot = (
                    row
                )
                mini[iz, it, im] = v_mini
                mact[iz, it, im] = v_mact
                logl[iz, it, im] = v_logl
                logt[iz, it, im] = v_logt
                logg[iz, it, im] = v_logg
                ffco[iz, it, im] = v_ffco
                phase[iz, it, im] = v_phase
                lmdot[iz, it, im] = v_lmdot

    tracks_grp = h5.require_group(f"/libraries/isochrones/{isoc_lib}/tracks")

    d_nmass = tracks_grp.create_dataset("nmass", data=nmass, dtype=np.int32)
    d_nmass.attrs["role"] = "isoc_nmass"
    d_nmass.attrs["dims_csv"] = "nt,nz"
    d_nmass.attrs["representation"] = "dense_nd"

    d_tstep = tracks_grp.create_dataset(
        "timestep", data=timestep_logyr, dtype=np.float64
    )
    d_tstep.attrs["role"] = "isoc_timestep"
    d_tstep.attrs["dims_csv"] = "nt,nz"
    d_tstep.attrs["representation"] = "dense_nd"

    def _write_3d(name: str, arr: np.ndarray) -> None:
        ds = tracks_grp.create_dataset(name, data=arr, dtype=np.float64)
        ds.attrs["role"] = f"isoc_{name}"
        ds.attrs["dims_csv"] = "nm,nt,nz"
        ds.attrs["representation"] = "dense_nd"

    _write_3d("mini", mini)
    _write_3d("mact", mact)
    _write_3d("logl", logl)
    _write_3d("logt", logt)
    _write_3d("logg", logg)
    _write_3d("phase", phase)
    _write_3d("ffco", ffco)
    _write_3d("lmdot", lmdot)


def _read_nebular_file_lines(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Missing nebular file: {path}")
    lines = [ln.strip() for ln in path.read_text().splitlines()]
    return [ln for ln in lines if ln]


def _interp_linear_extrap(
    x_in: np.ndarray, y_in: np.ndarray, x_out: np.ndarray
) -> np.ndarray:
    if x_in.ndim != 1 or y_in.ndim != 1 or x_out.ndim != 1:
        raise ValueError("Interpolation inputs must be 1D arrays")
    if x_in.size < 2 or y_in.size != x_in.size:
        raise ValueError("Interpolation input sizes are invalid")

    idx = np.searchsorted(x_in, x_out, side="right") - 1
    idx = np.clip(idx, 0, x_in.size - 2)

    x0 = x_in[idx]
    x1 = x_in[idx + 1]
    y0 = y_in[idx]
    y1 = y_in[idx + 1]
    slope = (y1 - y0) / (x1 - x0)
    return y0 + slope * (x_out - x0)


def _read_nebular_continuum(
    sps_home: Path,
    isoc_lib: str,
    prefix: str,
    axis_lambda: np.ndarray,
    n_z: int,
    n_age: int,
    n_u: int,
) -> np.ndarray:
    path = sps_home / "data" / "nebular" / f"ZAU_{prefix}_{isoc_lib}.cont"
    lines = _read_nebular_file_lines(path)
    if len(lines) < 2:
        raise ValueError(f"Nebular continuum file is malformed: {path}")

    raw_lam = np.fromstring(lines[1], sep=" ", dtype=np.float64)
    if raw_lam.size <= 0:
        raise ValueError(f"Nebular continuum wavelength row is empty: {path}")

    records_spec: list[np.ndarray] = []

    idx = 2
    while idx < len(lines):
        if idx + 1 >= len(lines):
            raise ValueError(
                f"Incomplete nebular continuum record near line {idx + 1}: {path}"
            )

        meta = np.fromstring(lines[idx], sep=" ", dtype=np.float64)
        spec = np.fromstring(lines[idx + 1], sep=" ", dtype=np.float64)
        idx += 2

        if meta.size < 3:
            raise ValueError(
                f"Invalid nebular continuum metadata row near line {idx - 1}: {path}"
            )
        if spec.size != raw_lam.size:
            raise ValueError(
                f"Continuum spectrum length mismatch near line {idx}: got {spec.size}, expected {raw_lam.size}"
            )

        records_spec.append(spec)

    if not records_spec:
        raise ValueError(f"No continuum records parsed from: {path}")

    n_lam = int(axis_lambda.size)
    expected_records = n_z * n_age * n_u
    if len(records_spec) != expected_records:
        raise ValueError(
            f"Nebular continuum record count mismatch: got {len(records_spec)}, expected {expected_records}"
        )

    cont = np.full((n_lam, n_z, n_age, n_u), np.float32(-1.0e30), dtype=np.float32)
    floor = 1.0e-95

    for rec, raw_spec in enumerate(records_spec):
        iz = rec // (n_age * n_u)
        rem = rec % (n_age * n_u)
        ia = rem // n_u
        iu = rem % n_u
        interp = _interp_linear_extrap(raw_lam, np.log10(raw_spec + floor), axis_lambda)
        cont[:, iz, ia, iu] = interp.astype(np.float32)

    return cont


def _read_nebular_lines(
    sps_home: Path,
    isoc_lib: str,
    prefix: str,
    n_z: int | None = None,
    n_age: int | None = None,
    n_u: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    path = sps_home / "data" / "nebular" / f"ZAU_{prefix}_{isoc_lib}.lines"
    lines = _read_nebular_file_lines(path)
    if len(lines) < 2:
        raise ValueError(f"Nebular lines file is malformed: {path}")

    line_pos = np.fromstring(lines[1], sep=" ", dtype=np.float64)
    if line_pos.size <= 0:
        raise ValueError(f"Nebular lines wavelength row is empty: {path}")

    metas: list[np.ndarray] = []
    spectra: list[np.ndarray] = []

    idx = 2
    floor = 1.0e-95
    while idx < len(lines):
        if idx + 1 >= len(lines):
            raise ValueError(
                f"Incomplete nebular lines record near line {idx + 1}: {path}"
            )

        meta = np.fromstring(lines[idx], sep=" ", dtype=np.float64)
        vals = np.fromstring(lines[idx + 1], sep=" ", dtype=np.float64)
        idx += 2

        if meta.size < 3:
            raise ValueError(
                f"Invalid nebular lines metadata row near line {idx - 1}: {path}"
            )
        if vals.size != line_pos.size:
            raise ValueError(
                f"Line spectrum length mismatch near line {idx}: got {vals.size}, expected {line_pos.size}"
            )

        metas.append(meta[:3])
        spectra.append(vals)

    if not metas:
        raise ValueError(f"No nebular lines records parsed from: {path}")

    meta_arr = np.asarray(metas, dtype=np.float64)
    n_records = meta_arr.shape[0]

    if n_z is None:
        n_z = int(np.unique(meta_arr[:, 0]).size)
    if n_age is None:
        n_age = int(np.unique(meta_arr[:, 1]).size)
    if n_u is None:
        n_u = int(np.unique(meta_arr[:, 2]).size)

    if n_records != n_z * n_age * n_u:
        raise ValueError(
            f"Nebular lines record count mismatch: got {n_records}, expected {n_z * n_age * n_u}"
        )

    line = np.full(
        (line_pos.size, n_z, n_age, n_u), np.float32(-1.0e30), dtype=np.float32
    )
    logz = np.zeros(n_z, dtype=np.float64)
    age = np.zeros(n_age, dtype=np.float64)
    logu = np.zeros(n_u, dtype=np.float64)

    for rec, (meta, vals) in enumerate(zip(meta_arr, spectra)):
        iz = rec // (n_age * n_u)
        rem = rec % (n_age * n_u)
        ia = rem // n_u
        iu = rem % n_u

        logz[iz] = meta[0]
        age[ia] = meta[1]
        logu[iu] = meta[2]
        line[:, iz, ia, iu] = np.log10(vals + floor).astype(np.float32)

    return logz, np.log10(age), logu, line_pos, line


def _write_nebular(
    h5: h5py.File, sps_home: Path, isoc_lib: str, axis_lambda: np.ndarray
) -> None:
    for prefix in ("WD", "ND"):
        logz, age_log, logu, line_pos, line = _read_nebular_lines(
            sps_home, isoc_lib, prefix
        )
        n_z = int(logz.size)
        n_age = int(age_log.size)
        n_u = int(logu.size)

        cont = _read_nebular_continuum(
            sps_home, isoc_lib, prefix, axis_lambda, n_z, n_age, n_u
        )

        grp = h5.require_group(f"/libraries/nebular/{prefix}")
        role_prefix = f"nebular_{prefix.lower()}"

        d_logz = grp.create_dataset("logz", data=logz, dtype=np.float64)
        d_logz.attrs["role"] = f"{role_prefix}_logz"
        d_logz.attrs["dims_csv"] = "z"
        d_logz.attrs["representation"] = "dense_nd"

        d_age = grp.create_dataset("age", data=age_log, dtype=np.float64)
        d_age.attrs["role"] = f"{role_prefix}_age"
        d_age.attrs["dims_csv"] = "age"
        d_age.attrs["representation"] = "dense_nd"

        d_logu = grp.create_dataset("logu", data=logu, dtype=np.float64)
        d_logu.attrs["role"] = f"{role_prefix}_logu"
        d_logu.attrs["dims_csv"] = "u"
        d_logu.attrs["representation"] = "dense_nd"

        d_line_pos = grp.create_dataset("line_pos", data=line_pos, dtype=np.float64)
        d_line_pos.attrs["role"] = f"{role_prefix}_line_pos"
        d_line_pos.attrs["dims_csv"] = "line"
        d_line_pos.attrs["representation"] = "dense_nd"

        # For Fortran rank-4 reads expecting (lam,z,age,u), write in C-order
        # as (u,age,z,lam).
        cont_c = np.transpose(cont, (3, 2, 1, 0)).astype(np.float32, copy=False)
        d_cont = grp.create_dataset("cont", data=cont_c, dtype=np.float32)
        d_cont.attrs["role"] = f"{role_prefix}_cont"
        d_cont.attrs["dims_csv"] = "lam,z,age,u"
        d_cont.attrs["representation"] = "dense_nd"

        # For Fortran rank-4 reads expecting (line,z,age,u), write in C-order
        # as (u,age,z,line).
        line_c = np.transpose(line, (3, 2, 1, 0)).astype(np.float32, copy=False)
        d_line = grp.create_dataset("lines", data=line_c, dtype=np.float32)
        d_line.attrs["role"] = f"{role_prefix}_line"
        d_line.attrs["dims_csv"] = "line,z,age,u"
        d_line.attrs["representation"] = "dense_nd"


def _read_wmbasic_file(path: Path, n_logt: int) -> tuple[np.ndarray, np.ndarray]:
    arr = np.loadtxt(path, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"Unexpected WMBasic array rank in {path}: {arr.ndim}")

    expected_cols = 1 + 3 * n_logt
    if arr.shape[1] != expected_cols:
        raise ValueError(
            f"Unexpected WMBasic column count in {path}: got {arr.shape[1]}, expected {expected_cols}"
        )

    lam = arr[:, 0]
    g1 = arr[:, 1 : 1 + n_logt]
    g2 = arr[:, 1 + n_logt : 1 + 2 * n_logt]
    g3 = arr[:, 1 + 2 * n_logt : 1 + 3 * n_logt]

    spec = np.stack((g1, g2, g3), axis=2)
    return lam, spec


def _read_wr_file(
    path: Path, n_logt: int, n_z: int = 5
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"WR file is empty: {path}")

    lam = np.fromstring(lines[0], sep=" ", dtype=np.float64)
    if lam.size <= 0:
        raise ValueError(f"WR wavelength row is empty: {path}")

    expected_payload = n_z * n_logt * 2
    payload = lines[1:]
    if len(payload) != expected_payload:
        raise ValueError(
            f"Unexpected WR payload line count in {path}: got {len(payload)}, expected {expected_payload}"
        )

    z = np.zeros(n_z, dtype=np.float64)
    spec = np.zeros((n_z, n_logt, lam.size), dtype=np.float32)

    k = 0
    for iz in range(n_z):
        for it in range(n_logt):
            meta = np.fromstring(payload[k], sep=" ", dtype=np.float64)
            k += 1
            if meta.size < 2:
                raise ValueError(
                    f"Invalid WR metadata row in {path} near payload line {k}"
                )
            if it == 0:
                z[iz] = meta[1]

            flux = np.fromstring(payload[k], sep=" ", dtype=np.float64)
            k += 1
            if flux.size != lam.size:
                raise ValueError(
                    f"Invalid WR spectrum row in {path}: got {flux.size} columns, expected {lam.size}"
                )
            spec[iz, it, :] = flux.astype(np.float32, copy=False)

    return lam, z, spec


def _write_auxiliary(h5: h5py.File, sps_home: Path) -> None:
    hot_dir = sps_home / "data" / "spectra" / "Hot_spectra"
    agb_dir = sps_home / "data" / "spectra" / "AGB_spectra"

    # ------------------------------
    # WMBasic
    # ------------------------------
    wmb_logt = np.loadtxt(hot_dir / "WMBASIC.teff", dtype=np.float64, ndmin=1)
    wmb_z = np.loadtxt(hot_dir / "WMBASIC_zlegend.dat", dtype=np.float64, ndmin=1)

    n_logt = int(wmb_logt.size)
    n_z = int(wmb_z.size)

    lam_ref: np.ndarray | None = None
    wmb_stack: list[np.ndarray] = []
    for z in wmb_z:
        zstr = f"{float(z):0.4f}"
        lam, spec = _read_wmbasic_file(hot_dir / f"WMBASIC_z{zstr}.spec", n_logt)
        if lam_ref is None:
            lam_ref = lam
        elif not np.allclose(lam, lam_ref, rtol=0.0, atol=0.0):
            raise ValueError("WMBasic wavelength grids differ between metallicities")
        wmb_stack.append(spec)

    if lam_ref is None:
        raise ValueError("No WMBasic spectra were parsed")

    # Shape: (n_z, n_lam, n_logt, n_logg)
    wmb_z_lam_t_g = np.stack(wmb_stack, axis=0)
    # Reorder to (n_z, n_logg, n_logt, n_lam) for Fortran read as (lam,logt,logg,z)
    wmb_c = np.transpose(wmb_z_lam_t_g, (0, 3, 2, 1)).astype(np.float32, copy=False)

    wmb_grp = h5.require_group("/libraries/auxiliary/wmbasic")

    ds = wmb_grp.create_dataset("logt", data=wmb_logt, dtype=np.float64)
    ds.attrs["role"] = "wmb_logt"
    ds.attrs["dims_csv"] = "logt"
    ds.attrs["representation"] = "dense_nd"

    ds = wmb_grp.create_dataset("z", data=wmb_z, dtype=np.float64)
    ds.attrs["role"] = "wmb_z"
    ds.attrs["dims_csv"] = "z"
    ds.attrs["representation"] = "dense_nd"

    ds = wmb_grp.create_dataset(
        "lam", data=lam_ref.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "wmb_lam"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = wmb_grp.create_dataset("spec", data=wmb_c, dtype=np.float32)
    ds.attrs["role"] = "wmb_spec"
    ds.attrs["dims_csv"] = "lam,logt,logg,z"
    ds.attrs["representation"] = "dense_nd"

    # ------------------------------
    # Post-AGB (Rauch)
    # ------------------------------
    pagb_teff = np.loadtxt(hot_dir / "ipagb.teff", dtype=np.float64, ndmin=1)
    pagb_logt = np.log10(pagb_teff)

    halo = np.loadtxt(hot_dir / "ipagb_halo.spec", dtype=np.float64)
    solar = np.loadtxt(hot_dir / "ipagb_solar.spec", dtype=np.float64)

    if halo.ndim != 2 or solar.ndim != 2:
        raise ValueError("Unexpected Post-AGB table rank")
    if halo.shape != solar.shape:
        raise ValueError("Post-AGB halo/solar table shapes differ")
    if halo.shape[1] != 1 + pagb_logt.size:
        raise ValueError(
            f"Unexpected Post-AGB column count: got {halo.shape[1]}, expected {1 + pagb_logt.size}"
        )

    pagb_lam = halo[:, 0]
    if not np.allclose(solar[:, 0], pagb_lam, rtol=0.0, atol=0.0):
        raise ValueError("Post-AGB halo/solar wavelength grids differ")

    # Raw shape: (n_lam, n_logt, n_z)
    pagb_raw = np.stack((halo[:, 1:], solar[:, 1:]), axis=2)
    # Reorder to (n_z, n_logt, n_lam) for Fortran read as (lam,logt,z)
    pagb_c = np.transpose(pagb_raw, (2, 1, 0)).astype(np.float32, copy=False)

    pagb_grp = h5.require_group("/libraries/auxiliary/pagb")

    ds = pagb_grp.create_dataset(
        "logt", data=pagb_logt.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "pagb_logt"
    ds.attrs["dims_csv"] = "logt"
    ds.attrs["representation"] = "dense_nd"

    ds = pagb_grp.create_dataset(
        "lam", data=pagb_lam.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "pagb_lam"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = pagb_grp.create_dataset("spec", data=pagb_c, dtype=np.float32)
    ds.attrs["role"] = "pagb_spec"
    ds.attrs["dims_csv"] = "lam,logt,z"
    ds.attrs["representation"] = "dense_nd"

    # ------------------------------
    # Wolf-Rayet (CMFGEN)
    # ------------------------------
    wrn_logt = np.loadtxt(hot_dir / "CMFGEN_WN.teff", dtype=np.float64, ndmin=1)
    wrc_logt = np.loadtxt(hot_dir / "CMFGEN_WC.teff", dtype=np.float64, ndmin=1)

    n_logt_wr = int(wrn_logt.size)
    if int(wrc_logt.size) != n_logt_wr:
        raise ValueError("WR WN/WC Teff grids have different sizes")

    wr_lam, wr_z, wrn_raw = _read_wr_file(
        hot_dir / "CMFGEN_WN_Zall.spec", n_logt_wr, n_z=5
    )
    wr_lam_wc, wr_z_wc, wrc_raw = _read_wr_file(
        hot_dir / "CMFGEN_WC_Zall.spec", n_logt_wr, n_z=5
    )

    if not np.allclose(wr_lam_wc, wr_lam, rtol=0.0, atol=0.0):
        raise ValueError("WR WN/WC wavelength grids differ")
    if not np.allclose(wr_z_wc, wr_z, rtol=0.0, atol=0.0):
        raise ValueError("WR WN/WC metallicity grids differ")

    # For Fortran rank-3 reads expecting (lam,logt,z), write C-order as (z,logt,lam).
    wrn_c = wrn_raw.astype(np.float32, copy=False)
    wrc_c = wrc_raw.astype(np.float32, copy=False)

    wr_grp = h5.require_group("/libraries/auxiliary/wr")

    ds = wr_grp.create_dataset(
        "logt_wn", data=wrn_logt.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "wr_logt_wn"
    ds.attrs["dims_csv"] = "logt"
    ds.attrs["representation"] = "dense_nd"

    ds = wr_grp.create_dataset(
        "logt_wc", data=wrc_logt.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "wr_logt_wc"
    ds.attrs["dims_csv"] = "logt"
    ds.attrs["representation"] = "dense_nd"

    ds = wr_grp.create_dataset("z", data=wr_z.astype(np.float64), dtype=np.float64)
    ds.attrs["role"] = "wr_z"
    ds.attrs["dims_csv"] = "z"
    ds.attrs["representation"] = "dense_nd"

    ds = wr_grp.create_dataset("lam", data=wr_lam.astype(np.float64), dtype=np.float64)
    ds.attrs["role"] = "wr_lam"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = wr_grp.create_dataset("spec_wn", data=wrn_c, dtype=np.float32)
    ds.attrs["role"] = "wr_spec_wn"
    ds.attrs["dims_csv"] = "lam,logt,z"
    ds.attrs["representation"] = "dense_nd"

    ds = wr_grp.create_dataset("spec_wc", data=wrc_c, dtype=np.float32)
    ds.attrs["role"] = "wr_spec_wc"
    ds.attrs["dims_csv"] = "lam,logt,z"
    ds.attrs["representation"] = "dense_nd"

    # ------------------------------
    # AGB (Lancon+Wood, Aringer)
    # ------------------------------
    orich_teff = np.loadtxt(agb_dir / "Orich.teff", dtype=np.float64, skiprows=1)
    if orich_teff.ndim != 2 or orich_teff.shape[1] != 23:
        raise ValueError("Unexpected O-rich Teff table shape")

    agb_z_o = orich_teff[0, 1:]
    agb_logt_o = np.log10(orich_teff[1:, 1:].T)

    crich_teff = np.loadtxt(agb_dir / "Crich.teff", dtype=np.float64, skiprows=1)
    if crich_teff.ndim != 2 or crich_teff.shape[1] != 2:
        raise ValueError("Unexpected C-rich Teff table shape")
    agb_logt_c = np.log10(crich_teff[:, 1])

    aringer_teff = np.loadtxt(agb_dir / "Crich_Aringer.teff", dtype=np.float64, ndmin=1)
    agb_logt_car = np.log10(aringer_teff)

    orich_spec = np.loadtxt(agb_dir / "Orich.spec", dtype=np.float64)
    crich_spec = np.loadtxt(agb_dir / "Crich.spec", dtype=np.float64)
    aringer_spec = np.loadtxt(agb_dir / "Crich_Aringer.spec", dtype=np.float64)

    if orich_spec.ndim != 2 or crich_spec.ndim != 2 or aringer_spec.ndim != 2:
        raise ValueError("Unexpected AGB spectra table rank")

    lam_o = orich_spec[:, 0]
    lam_c = crich_spec[:, 0]
    lam_car = aringer_spec[:, 0]

    spec_o = orich_spec[:, 1:]
    spec_c = crich_spec[:, 1:]
    spec_car = aringer_spec[:, 1:]

    if spec_o.shape[1] != agb_logt_o.shape[1]:
        raise ValueError("O-rich AGB spectra/Teff size mismatch")
    if spec_c.shape[1] != agb_logt_c.size:
        raise ValueError("C-rich AGB spectra/Teff size mismatch")
    if spec_car.shape[1] != agb_logt_car.size:
        raise ValueError("Aringer AGB spectra/Teff size mismatch")

    # For Fortran rank-2 reads expecting (z,logt), write C-order as (logt,z).
    logt_o_c = np.transpose(agb_logt_o, (1, 0)).astype(np.float64, copy=False)
    # For Fortran rank-2 reads expecting (lam,logt), write C-order as (logt,lam).
    spec_o_c = np.transpose(spec_o, (1, 0)).astype(np.float32, copy=False)
    spec_c_c = np.transpose(spec_c, (1, 0)).astype(np.float32, copy=False)
    spec_car_c = np.transpose(spec_car, (1, 0)).astype(np.float32, copy=False)

    agb_grp = h5.require_group("/libraries/auxiliary/agb")

    ds = agb_grp.create_dataset(
        "z_o", data=agb_z_o.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "agb_z_o"
    ds.attrs["dims_csv"] = "z"
    ds.attrs["representation"] = "dense_nd"

    ds = agb_grp.create_dataset(
        "logt_c", data=agb_logt_c.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "agb_logt_c"
    ds.attrs["dims_csv"] = "logt"
    ds.attrs["representation"] = "dense_nd"

    ds = agb_grp.create_dataset(
        "logt_car", data=agb_logt_car.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "agb_logt_car"
    ds.attrs["dims_csv"] = "logt"
    ds.attrs["representation"] = "dense_nd"

    ds = agb_grp.create_dataset(
        "lam_o", data=lam_o.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "agb_lam_o"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = agb_grp.create_dataset(
        "lam_c", data=lam_c.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "agb_lam_c"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = agb_grp.create_dataset(
        "lam_car", data=lam_car.astype(np.float64), dtype=np.float64
    )
    ds.attrs["role"] = "agb_lam_car"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = agb_grp.create_dataset("logt_o", data=logt_o_c, dtype=np.float64)
    ds.attrs["role"] = "agb_logt_o"
    ds.attrs["dims_csv"] = "z,logt"
    ds.attrs["representation"] = "dense_nd"

    ds = agb_grp.create_dataset("spec_o", data=spec_o_c, dtype=np.float32)
    ds.attrs["role"] = "agb_spec_o"
    ds.attrs["dims_csv"] = "lam,logt"
    ds.attrs["representation"] = "dense_nd"

    ds = agb_grp.create_dataset("spec_c", data=spec_c_c, dtype=np.float32)
    ds.attrs["role"] = "agb_spec_c"
    ds.attrs["dims_csv"] = "lam,logt"
    ds.attrs["representation"] = "dense_nd"

    ds = agb_grp.create_dataset("spec_car", data=spec_car_c, dtype=np.float32)
    ds.attrs["role"] = "agb_spec_car"
    ds.attrs["dims_csv"] = "lam,logt"
    ds.attrs["representation"] = "dense_nd"


def _normalize_dust_type(dust_type: str) -> str:
    return "THEMIS" if dust_type.strip().lower() == "themis" else "DL07"


def _write_dust_emission(h5: h5py.File, sps_home: Path, dust_type: str) -> None:
    dust = _normalize_dust_type(dust_type)
    dust_dir = sps_home / "data" / "dust" / "dustem"

    if dust == "THEMIS":
        qpah = QPAH_ARR_THEMIS
        umin = UMIN_ARR_THEMIS
        n_lam = 576
    else:
        qpah = QPAH_ARR_DL07
        umin = UMIN_ARR_DL07
        n_lam = 1001

    n_qpah = int(qpah.size)
    n_umin_cols = int(umin.size) * 2

    lam_ref: np.ndarray | None = None
    # Raw target shape for Fortran (lam, qpah, umin_cols): write C-order as (umin_cols, qpah, lam)
    spec_c = np.zeros((n_umin_cols, n_qpah, n_lam), dtype=np.float32)

    for k in range(n_qpah):
        if k == 10:
            path = dust_dir / f"{dust}_MW3.1_100.dat"
        else:
            path = dust_dir / f"{dust}_MW3.1_{k}0.dat"

        table = np.loadtxt(path, dtype=np.float64, skiprows=2)
        if table.ndim != 2 or table.shape[1] != 1 + n_umin_cols:
            raise ValueError(
                f"Unexpected dust emission table shape in {path}: got {table.shape}, expected (*, {1 + n_umin_cols})"
            )

        lam = table[:, 0] * 1.0e4  # microns -> Angstroms
        if lam.size != n_lam:
            raise ValueError(
                f"Unexpected dust emission wavelength count in {path}: got {lam.size}, expected {n_lam}"
            )
        if lam_ref is None:
            lam_ref = lam
        elif not np.allclose(lam, lam_ref, rtol=0.0, atol=0.0):
            raise ValueError("Dust emission wavelength grids differ between qpah files")

        # table[:, 1:] is (lam, umin_cols) -> transpose to (umin_cols, lam)
        spec_c[:, k, :] = table[:, 1:].T.astype(np.float32, copy=False)

    if lam_ref is None:
        raise ValueError("No dust emission templates were parsed")

    grp = h5.require_group(f"/libraries/dust/emission/{dust}")

    ds = grp.create_dataset("qpah", data=qpah.astype(np.float64), dtype=np.float64)
    ds.attrs["role"] = "dust_em_qpah"
    ds.attrs["dims_csv"] = "qpah"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("umin", data=umin.astype(np.float64), dtype=np.float64)
    ds.attrs["role"] = "dust_em_umin"
    ds.attrs["dims_csv"] = "umin"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("lam", data=lam_ref.astype(np.float64), dtype=np.float64)
    ds.attrs["role"] = "dust_em_lam"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("spec", data=spec_c, dtype=np.float32)
    ds.attrs["role"] = "dust_em_spec"
    ds.attrs["dims_csv"] = "lam,qpah,umin"
    ds.attrs["representation"] = "dense_nd"


def _write_agn_dust(h5: h5py.File, sps_home: Path) -> None:
    path = sps_home / "data" / "dust" / "Nenkova08_y010_torusg_n10_q2.0.dat"

    with path.open("r", encoding="utf-8") as f:
        # header
        for _ in range(3):
            _ = f.readline()
        tau = np.fromstring(f.readline(), sep=" ", dtype=np.float64)
        rows = [
            np.fromstring(line, sep=" ", dtype=np.float64) for line in f if line.strip()
        ]

    if tau.size == 0:
        raise ValueError(f"Failed to parse AGN tau grid from {path}")
    if not rows:
        raise ValueError(f"Failed to parse AGN spectra table from {path}")

    table = np.vstack(rows)
    if table.ndim != 2 or table.shape[1] != 1 + tau.size:
        raise ValueError(
            f"Unexpected AGN dust table shape in {path}: got {table.shape}, expected (*, {1 + tau.size})"
        )

    lam = table[:, 0].astype(np.float64, copy=False)
    # Raw (lam, tau) -> C-order for Fortran (lam, tau) is (tau, lam)
    spec_c = table[:, 1:].T.astype(np.float32, copy=False)

    grp = h5.require_group("/libraries/dust/agn")

    ds = grp.create_dataset("tau", data=tau.astype(np.float64), dtype=np.float64)
    ds.attrs["role"] = "agn_dust_tau"
    ds.attrs["dims_csv"] = "tau"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("lam", data=lam, dtype=np.float64)
    ds.attrs["role"] = "agn_dust_lam"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("spec", data=spec_c, dtype=np.float32)
    ds.attrs["role"] = "agn_dust_spec"
    ds.attrs["dims_csv"] = "lam,tau"
    ds.attrs["representation"] = "dense_nd"


def _write_attenuation(h5: h5py.File, sps_home: Path) -> None:
    dust_dir = sps_home / "data" / "dust"

    n_lam = 25
    n_geom = 18
    n_tau = 6

    shell = np.loadtxt(dust_dir / "alldirty_h.dat", dtype=np.float64, skiprows=3)
    cloudy = np.loadtxt(dust_dir / "alldirty_c.dat", dtype=np.float64, skiprows=3)

    expected_rows = n_geom * n_lam
    expected_cols = 2 + n_tau
    if shell.ndim != 2 or shell.shape != (expected_rows, expected_cols):
        raise ValueError(
            f"Unexpected WG shell table shape: got {shell.shape}, expected {(expected_rows, expected_cols)}"
        )
    if cloudy.ndim != 2 or cloudy.shape != (expected_rows, expected_cols):
        raise ValueError(
            f"Unexpected WG cloudy table shape: got {cloudy.shape}, expected {(expected_rows, expected_cols)}"
        )

    wg_lam = shell[:n_lam, 0].astype(np.float64, copy=False)
    wg_spec_c = np.zeros((2, n_tau, n_geom, n_lam), dtype=np.float32)

    for geom in range(n_geom):
        sl = slice(geom * n_lam, (geom + 1) * n_lam)
        lam_s = shell[sl, 0]
        lam_c = cloudy[sl, 0]
        if not np.allclose(lam_s, wg_lam, rtol=0.0, atol=1.0e-12):
            raise ValueError("WG shell wavelength grid mismatch between geometries")
        if not np.allclose(lam_c, wg_lam, rtol=0.0, atol=1.0e-12):
            raise ValueError("WG cloudy wavelength grid mismatch")

        wg_spec_c[0, :, geom, :] = shell[sl, 2:].T.astype(np.float32, copy=False)
        wg_spec_c[1, :, geom, :] = cloudy[sl, 2:].T.astype(np.float32, copy=False)

    smc = np.loadtxt(dust_dir / "Gordon03_table4.dat", dtype=np.float64)
    if smc.ndim != 2 or smc.shape[1] < 3:
        raise ValueError("Unexpected SMC attenuation table shape")

    smc_lam = (smc[:, 0][::-1] * 1.0e4).astype(np.float64, copy=False)
    smc_ext = smc[:, 2][::-1].astype(np.float64, copy=False)

    grp = h5.require_group("/libraries/dust/attenuation")

    ds = grp.create_dataset("wg_lam", data=wg_lam, dtype=np.float64)
    ds.attrs["role"] = "dust_att_wg_lam"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("wg_spec", data=wg_spec_c, dtype=np.float32)
    ds.attrs["role"] = "dust_att_wg_spec"
    ds.attrs["dims_csv"] = "lam,geom,tau,type"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("smc_lam", data=smc_lam, dtype=np.float64)
    ds.attrs["role"] = "dust_att_smc_lam"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("smc_ext", data=smc_ext, dtype=np.float64)
    ds.attrs["role"] = "dust_att_smc_ext"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"


def _write_xrb(h5: h5py.File, sps_home: Path) -> None:
    xrb_dir = sps_home / "data" / "spectra" / "xrb"

    lam = np.loadtxt(xrb_dir / "xsp.lambda", dtype=np.float64, ndmin=1)
    age = (
        np.log10(
            np.array(
                [1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0, 12.6, 16.0, 20.0], dtype=np.float64
            )
        )
        + 6.0
    )
    z = np.array(
        [-1.3, -1.0, -0.8, -0.7, -0.5, -0.4, -0.3, -0.2, 0.0, 0.2, 0.3],
        dtype=np.float64,
    )
    z_tags = [
        "-1.30",
        "-1.00",
        "-0.80",
        "-0.70",
        "-0.50",
        "-0.40",
        "-0.30",
        "-0.20",
        "+0.00",
        "+0.20",
        "+0.30",
    ]

    n_lam = int(lam.size)
    n_age = int(age.size)
    n_z = int(z.size)

    spec_c = np.zeros((n_z, n_age, n_lam), dtype=np.float32)
    for iz, tag in enumerate(z_tags):
        arr = np.loadtxt(xrb_dir / f"xsp_feh{tag}.spec", dtype=np.float64)
        if arr.ndim == 1:
            arr = arr[np.newaxis, :]
        if arr.shape != (n_age, n_lam):
            raise ValueError(
                f"Unexpected XRB spectra shape for [Fe/H]={tag}: got {arr.shape}, expected {(n_age, n_lam)}"
            )
        spec_c[iz, :, :] = arr.astype(np.float32, copy=False)

    grp = h5.require_group("/libraries/xrb")

    ds = grp.create_dataset("lam", data=lam.astype(np.float64), dtype=np.float64)
    ds.attrs["role"] = "xrb_lam"
    ds.attrs["dims_csv"] = "lam"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("age", data=age.astype(np.float64), dtype=np.float64)
    ds.attrs["role"] = "xrb_age"
    ds.attrs["dims_csv"] = "age"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("z", data=z.astype(np.float64), dtype=np.float64)
    ds.attrs["role"] = "xrb_z"
    ds.attrs["dims_csv"] = "z"
    ds.attrs["representation"] = "dense_nd"

    ds = grp.create_dataset("spec", data=spec_c, dtype=np.float32)
    ds.attrs["role"] = "xrb_spec"
    ds.attrs["dims_csv"] = "lam,age,z"
    ds.attrs["representation"] = "dense_nd"


def convert(
    sps_home: Path, spec_lib: str, isoc_lib: str, out_file: Path, dust_type: str
) -> None:
    spectra_dir = _spectra_dir(sps_home, spec_lib)
    spec = spec_lib.lower()

    if spec == "bpass":
        axis_lambda = _read_axis(spectra_dir / "bpass.lambda")
        axis_z = _read_axis(spectra_dir / "zlegend.dat")
        axis_logt, _mass_ssp = _read_bpass_mass_table(spectra_dir / "bpass.mass", int(axis_z.size))
        axis_logg = np.array([0.0], dtype=np.float64)

        cube_l_t_z = _read_bpass_cube(
            spectra_dir / "bpass_v2.2_salpeter100.ssp.bin",
            int(axis_lambda.size),
            int(axis_logt.size),
            int(axis_z.size),
        )
        cube_l_z_t_g = np.transpose(cube_l_t_z, (0, 2, 1))[:, :, :, np.newaxis]
        spectral_grid_c = np.transpose(cube_l_z_t_g, (3, 2, 1, 0)).astype(
            np.float32, copy=False
        )
    else:
        axis_lambda = _read_axis(spectra_dir / _lambda_filename(spec_lib))
        axis_z = _read_axis(spectra_dir / "zlegend.dat")

        # FSPS legacy spectral grid axes are shared from BaSeL tables.
        basel_dir = sps_home / "data" / "spectra" / "BaSeL3.1"
        axis_logt = _read_axis(basel_dir / "basel_logt.dat")
        axis_logg = _read_axis(basel_dir / "basel_logg.dat")

        n_lam = axis_lambda.size
        n_z = axis_z.size
        n_logt = axis_logt.size
        n_logg = axis_logg.size

        cubes = []
        for iz in range(n_z):
            zval = float(axis_z[iz])
            bin_path = _binary_filename(spec_lib, zval, spectra_dir)
            cube = _read_legacy_cube(bin_path, n_lam, n_logt, n_logg)
            cubes.append(cube)

        # Shape after stacking: (n_z, n_logg, n_logt, n_lam)
        stacked = np.stack(cubes, axis=0)

        # For Fortran backend expecting flux(lambda, z, logt, logg), write the dataset
        # in C-order as (n_logg, n_logt, n_z, n_lambda).
        spectral_grid_c = np.transpose(stacked, (1, 2, 0, 3)).astype(np.float32, copy=False)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_file, "w") as h5:
        axes_grp = h5.require_group("/axes")
        axes_grp.create_dataset("lambda", data=axis_lambda.astype(np.float64))
        axes_grp.create_dataset("z", data=axis_z.astype(np.float64))
        axes_grp.create_dataset("logt", data=axis_logt.astype(np.float64))
        axes_grp.create_dataset("logg", data=axis_logg.astype(np.float64))

        base_grp = h5.require_group(f"/libraries/spectra/{spec_lib}/base")
        dset = base_grp.create_dataset(
            "spectral_grid_nd", data=spectral_grid_c, dtype=np.float32
        )

        dset.attrs["role"] = "spectral_base"
        dset.attrs["dims_csv"] = "lambda,z,logt,logg"
        dset.attrs["unit"] = "Lsun/Hz/Msun"
        dset.attrs["has_missing_value"] = np.int32(0)

        _write_isochrones(h5, sps_home, isoc_lib)
        _write_nebular(h5, sps_home, isoc_lib, axis_lambda)
        _write_auxiliary(h5, sps_home)
        _write_dust_emission(h5, sps_home, dust_type)
        _write_agn_dust(h5, sps_home)
        _write_attenuation(h5, sps_home)
        _write_xrb(h5, sps_home)


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert legacy FSPS spectral binaries to FSDS HDF5."
    )
    p.add_argument("--sps-home", required=True, help="FSPS root path.")
    p.add_argument(
        "--spec-lib",
        required=True,
        help="Spectral library name, e.g. miles, basel, bpass, c3k_afe+0.0.",
    )
    p.add_argument(
        "--isoc-lib",
        required=True,
        help="Isochrone library name, e.g. mist, pdva, prsc, bsti, gnva, bpss.",
    )
    p.add_argument(
        "--dust-type",
        default="DL07",
        help="Dust emission library name: DL07 or THEMIS.",
    )
    p.add_argument("--out", required=True, help="Output HDF5 filename.")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    convert(
        Path(args.sps_home),
        args.spec_lib,
        args.isoc_lib,
        Path(args.out),
        args.dust_type,
    )
