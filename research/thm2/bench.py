#!/usr/bin/env python3
"""Experiment driver: python bench.py <mode>   (modes: check, scale_n, scale_m, base)"""
import sys, math, json, time
import thm2
from thm2 import run, fx


def line(r):
    st = r["st"]
    s = (f"n={r['n']:>7} m_bits={r['m_bits']:>7} M={r['M']:>5} N={r['N']:>6} "
         f"err={r.get('err', float('nan')):.1e} tB={r['tB']:7.2f} tC={r['tC']:8.2f} "
         f"art={st.t.get('art', 0):7.2f} luc={st.t.get('lucas', 0):5.2f} swp={st.t.get('sweep', 0):5.2f} "
         f"chunks={st.chunks:>5} leaf={st.leaf_steps:>11} sweep={st.sweep_steps:>8} "
         f"peak~{st.peak_bits}")
    if "tC1" in r:
        s += f" | thm1 tC={r['tC1']:.2f}s ops={r['ops1']} agree={r['agree']}"
    return s


def main():
    mode = sys.argv[1]
    rows = []
    if mode == "check":
        for n, mb in [(300, 128), (1000, 300), (3000, 700), (5000, 1000), (10000, 2000),
                      (10000, 500), (20000, 3000), (50000, 4000), (100000, 8000)]:
            r = run(n, mb, baseline=(n <= 5000))
            print(line(r), flush=True)
    elif mode == "base":
        for n in [1000, 2000, 4000, 8000, 16000]:
            mb = int(4 * math.sqrt(n) * math.log2(10))
            r = run(n, mb, baseline=True, check=False)
            print(line(r), flush=True)
    elif mode == "scale_n":
        # m ~ sqrt(n) decimal digits (x const), the headline case of Theorem 2
        c = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
        for n in [10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000]:
            mb = int(c * math.sqrt(n) * math.log2(10))
            r = run(n, mb, check=(n <= 1280000))
            print(line(r), flush=True)
    elif mode == "scale_m":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 40000
        for mb in [512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]:
            r = run(n, mb, check=False)
            print(line(r), flush=True)
    elif mode == "fixedm":
        # fixed memory, growing n: should look ~quadratic (up to logs)
        mb = int(sys.argv[2])
        for n in [10000, 20000, 40000, 80000, 160000, 320000]:
            r = run(n, mb, check=False)
            print(line(r), flush=True)
    elif mode == "thm1c":
        # complete Theorem-1 algorithm in C at Gourdon's own parameters
        import subprocess
        exe = sys.argv[2]
        for n in [int(x) for x in sys.argv[3].split(",")]:
            M = 2 * math.ceil(n / math.log(n) ** 3)
            M = max(M, 4)
            N = thm2.choose_N(n, 12, M)
            out = subprocess.run([exe, str(n), str(M), str(N)], capture_output=True, text=True).stdout.strip()
            f = float(out.split()[0].split("=")[1])
            ref = thm2.reference(n) / thm2.ONE
            print(f"n={n:>8} M={M:>5} N={N:>7} {out} err={abs(f-ref):.1e}", flush=True)
    elif mode == "fixedM":
        # hold M (hence N) fixed, vary memory: isolates the C-part dependence on m
        n = int(sys.argv[2]); M = int(sys.argv[3])
        for mb in [256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072]:
            r = run(n, mb, M=M, check=False)
            print(line(r), flush=True)


if __name__ == "__main__":
    main()
