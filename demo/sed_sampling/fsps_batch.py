"""Minimal ctypes wrapper of libfsps's context C API (include/fsps.h) for the batched CSP.

Only what the demo needs: a context, its settings, its grids, and
fsps_context_compute_csp_batch. All arrays cross the boundary as float64; the batch's outputs
are column-major (quantity x n) in C, i.e. C-ordered numpy arrays of shape (n, quantity).
"""

import ctypes as ct
import os

import numpy as np

_dbl_p = ct.POINTER(ct.c_double)
_int_p = ct.POINTER(ct.c_int)


def _ptr(a):
    return None if a is None else a.ctypes.data_as(_dbl_p)


class FspsError(RuntimeError):
    pass


class Context:
    """One FSPS context (its own settings, libraries and device copy)."""

    def __init__(self, lib_path, isoc="mist", spec="miles", dust="DL07", ints=None, floats=None):
        """ints / floats: settings applied before fsps_context_setup (and again after it). The
        cosmology must be set before setup: the redshift / distance table is built there, and
        H0, om0, ol0 set afterwards do not reach it."""
        self.lib = lib = ct.CDLL(os.path.abspath(lib_path))
        lib.fsps_context_create.argtypes = [_int_p, _int_p]
        lib.fsps_context_destroy.argtypes = [ct.c_int, _int_p]
        lib.fsps_context_setup.argtypes = [ct.c_int, ct.c_char_p, ct.c_char_p, ct.c_char_p, ct.c_int, _int_p]
        lib.fsps_context_set_int.argtypes = [ct.c_int, ct.c_char_p, ct.c_int, _int_p]
        lib.fsps_context_set_float.argtypes = [ct.c_int, ct.c_char_p, ct.c_double, _int_p]
        for name in ("nspec", "ntfull", "nbands", "nindx", "nz", "nemline"):
            getattr(lib, "fsps_context_get_" + name).argtypes = [ct.c_int, _int_p, _int_p]
        for name in ("timefull", "lambda", "filter_leff"):
            getattr(lib, "fsps_context_get_" + name).argtypes = [ct.c_int, _dbl_p, _int_p]
        lib.fsps_get_last_error.argtypes = [_int_p, ct.c_char_p, ct.c_int]
        lib.fsps_context_compute_csp_batch.argtypes = [
            ct.c_int, ct.c_int, ct.c_int, ct.POINTER(ct.c_char_p), ct.POINTER(_dbl_p), _dbl_p,
            _dbl_p, _dbl_p, _dbl_p, _dbl_p, _dbl_p, _dbl_p, _dbl_p, _dbl_p, _int_p]

        self._st = ct.c_int(0)
        h = ct.c_int(0)
        lib.fsps_context_create(ct.byref(h), ct.byref(self._st))
        self._check("create")
        self.handle = h.value
        for k, v in (ints or {}).items():
            self.set_int(k, v)
        for k, v in (floats or {}).items():
            self.set_float(k, v)
        lib.fsps_context_setup(-1, isoc.encode(), spec.encode(), dust.encode(), self.handle, ct.byref(self._st))
        self._check("setup")
        for k, v in (ints or {}).items():
            self.set_int(k, v)
        for k, v in (floats or {}).items():
            self.set_float(k, v)

    def _check(self, what):
        if self._st.value != 0:
            msg = ct.create_string_buffer(256)
            st = ct.c_int(0)
            self.lib.fsps_get_last_error(ct.byref(st), msg, 256)
            raise FspsError(f"{what}: status {self._st.value}: {msg.value.decode(errors='replace').strip()}")

    def set_int(self, key, value):
        self.lib.fsps_context_set_int(self.handle, key.encode(), int(value), ct.byref(self._st))
        self._check(f"set_int({key})")

    def set_float(self, key, value):
        self.lib.fsps_context_set_float(self.handle, key.encode(), float(value), ct.byref(self._st))
        self._check(f"set_float({key})")

    def dim(self, name):
        v = ct.c_int(0)
        getattr(self.lib, "fsps_context_get_" + name)(self.handle, ct.byref(v), ct.byref(self._st))
        self._check(name)
        return v.value

    def _grid(self, name, n):
        a = np.empty(n)
        getattr(self.lib, "fsps_context_get_" + name)(self.handle, _ptr(a), ct.byref(self._st))
        self._check(name)
        return a

    def log_age_grid(self):
        """SSP ages, log10(yr)."""
        return self._grid("timefull", self.dim("ntfull"))

    def wavelengths(self):
        return self._grid("lambda", self.dim("nspec"))

    def filter_leff(self):
        return self._grid("filter_leff", self.dim("nbands"))

    def compute_csp_batch(self, params, ssp_weights, mags, spec=None, mass=None, lbol=None, mdust=None,
                          emlines=None, indx=None):
        """One batched call: params maps a parameter name to its n per-vector values; ssp_weights
        (n, ntfull) is each vector's mass formed per SSP age; outputs are written in place."""
        n = mags.shape[0]
        names = list(params)
        vals = [np.ascontiguousarray(params[k], dtype=np.float64) for k in names]
        for v in vals:
            if v.shape != (n,):
                raise ValueError("every parameter array must have length n")
        c_names = (ct.c_char_p * len(names))(*[k.encode() for k in names])
        c_vals = (_dbl_p * len(names))(*[_ptr(v) for v in vals])
        self.lib.fsps_context_compute_csp_batch(
            self.handle, n, len(names), c_names, c_vals, _ptr(ssp_weights), _ptr(mags), _ptr(spec),
            _ptr(mass), _ptr(lbol), _ptr(mdust), None, _ptr(emlines), _ptr(indx), ct.byref(self._st))
        self._check("compute_csp_batch")

    def close(self):
        if self.handle:
            self.lib.fsps_context_destroy(self.handle, ct.byref(self._st))
            self.handle = 0
