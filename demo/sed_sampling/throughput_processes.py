#!/usr/bin/env python3
"""Aggregate throughput of P processes sharing the GPU, each running demo_batch_sed's model on
its own fresh draws (as P sampler workers would). One process alternates host work (SFH
weights, parameter sets, set-up of each call) with GPU work; a second one fills the GPU while
the first is on the host. The figure is total wall time / total vectors, measured from the
first process's start to the last one's end, after warm-up.

Usage: throughput_processes.py LIB [N [P [CALLS]]]   (defaults 10000, 2, 10)
"""

import multiprocessing as mp
import sys
import time

import numpy as np


def worker(lib, n, calls, seed, barrier, q):
    import demo_batch_sed as d
    model = d.BatchModel(lib)
    rng = np.random.default_rng(seed)
    out = model.outputs(n)
    for _ in range(2):
        model.evaluate(d.draw_vectors(rng, n), out)   # warm-up (SSP bases, device copy)
    proposals = [d.draw_vectors(rng, n) for _ in range(calls)]
    barrier.wait()
    t0 = time.perf_counter()
    for p in proposals:
        model.evaluate(p, out)
    q.put((t0, time.perf_counter()))


def main():
    lib = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    calls = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    ctx = mp.get_context("spawn")
    barrier, q = ctx.Barrier(nproc), ctx.Queue()
    procs = [ctx.Process(target=worker, args=(lib, n, calls, s, barrier, q)) for s in range(nproc)]
    for p in procs:
        p.start()
    spans = [q.get() for _ in procs]
    for p in procs:
        p.join()
    wall = max(t1 for _, t1 in spans) - min(t0 for t0, _ in spans)
    total = nproc * calls * n
    print(f"{nproc} process(es), N = {n} per call: {1e6 * wall / total:.2f} us per vector aggregate "
          f"({total} vectors in {wall:.2f} s)")


if __name__ == "__main__":
    main()
