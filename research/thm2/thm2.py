#!/usr/bin/env python3
"""Reconstruction of Gourdon's Theorem 2 (n-th decimal digit of pi with O(m) memory).

frac(10^n pi) ~= frac(B - C) with (Gourdon 2003, Prop. 1)
    B = sum_{k<(M+1)N} (-1)^k (4*10^n mod (2k+1)) / (2k+1)
    C = sum_{k<N}      (-1)^k (X s_k mod m_k) / m_k,
    X = 5^(N-2) 10^(n-N+2),  m_k = 2MN+2k+1,  s_k = sum_{j<=k} binom(N, j).

Theorem 1 computes every s_k mod m_k separately in O(k) word operations.
Here the C part is done with a *chunked accumulating remainder tree*:
all needed (target t, modulus q) pairs are sorted by t and cut into chunks
whose modulus product Q has ~m bits; for each chunk the prefix product of the
2x2 recurrence matrices up to the chunk's first target is formed by binary
splitting (exact products of ~m bits, then reduced mod Q), and the chunk's own
targets are resolved by a remainder-tree descent.

The division problem (prime factors p <= t of m_k are not invertible against
t!) is removed by splitting each m_k into coprime prime-power parts and using
partial fractions, frac(X s / m) = sum_q frac(X s ((m/q)^-1 mod q) / q):
  * part with primes p > t          -> ART item (target t, modulus q)       [main]
  * p <= t, p || m_k                -> Lucas: needs s_r, C(N,r) mod p with
                                       r = t mod p < p -> ART item (r, p)    [lucas]
  * p <= t, p^e || m_k, e >= 2      -> p-adic Lucas recursion (PadicBinom),
                                       O(e^2 p log_p N) word ops per target  [padic, default]
    (alternative small='sweep': p <= P_s one O(N) sweep per prime; p > P_s
     ART item (t, p^(e+v_p(t!))) -- kept for cross-checking)
Symmetry s_k = 2^N - s_{N-1-k} keeps every target t <= N/2.
"""
import math, sys, time
from gmpy2 import mpz, invert, f_mod
import numpy as np

FRAC_BITS = 96
ONE = 1 << FRAC_BITS


# ----------------------------------------------------------------- parameters
def choose_N(n, n0, M):
    N = math.ceil((n + n0 + 1) * math.log(10) / math.log(2 * math.e * M))
    N += N & 1
    assert N <= n + 2, "need n >= ~4 n0"
    return N


# ----------------------------------------------------------------- B part
def powmod_vec(base, e, mods):
    """base^e mod mods, vectorised; mods < 2^32."""
    mods = mods.astype(np.uint64)
    r = np.ones_like(mods)
    b = np.uint64(base) % mods
    while e:
        if e & 1:
            r = (r * b) % mods
        b = (b * b) % mods
        e >>= 1
    return r


def b_part(n, M, N, block=1 << 20):
    """frac(B) as a FRAC_BITS fixed-point integer (float64 terms, fsum)."""
    K = (M + 1) * N
    assert 2 * K + 1 < 1 << 32
    tot = 0.0
    parts = []
    for lo in range(0, K, block):
        k = np.arange(lo, min(K, lo + block), dtype=np.uint64)
        q = 2 * k + 1
        x = (4 * powmod_vec(10, n, q)) % q
        v = x.astype(np.float64) / q.astype(np.float64)
        v[(k & np.uint64(1)) == 1] *= -1
        parts.append(math.fsum(v))
    tot = math.fsum(parts)
    tot -= math.floor(tot)
    return int(tot * ONE) % ONE


# ----------------------------------------------------------------- number theory helpers
def primes_upto(x):
    s = bytearray([1]) * (x + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(x ** 0.5) + 1):
        if s[i]:
            s[i * i :: i] = bytearray(len(range(i * i, x + 1, i)))
    return [i for i in range(x + 1) if s[i]]


def factor_interval(M, N):
    """Factorisations of m_k = 2MN+2k+1, k<N (segmented sieve).
    NOTE: prototype keeps all N factorisations; a memory-faithful version
    sieves chunk by chunk (see docs)."""
    base = 2 * M * N + 1
    top = base + 2 * (N - 1)
    rem = [base + 2 * k for k in range(N)]
    fac = [[] for _ in range(N)]
    for p in primes_upto(math.isqrt(top) + 1):
        if p == 2:
            continue
        # first k with p | base + 2k
        k0 = (-base * pow(2, -1, p)) % p
        for k in range(k0, N, p):
            e = 0
            while rem[k] % p == 0:
                rem[k] //= p
                e += 1
            fac[k].append((p, e))
    for k in range(N):
        if rem[k] > 1:
            fac[k].append((rem[k], 1))
    return fac


def vp_fact(t, p):
    v, q = 0, p
    while q <= t:
        v += t // q
        q *= p
    return v


# ----------------------------------------------------------------- binary splitting of the recurrence
# state at target t:  P_t = N!/(N-t)!,  D_t = t!,  T_t = t! * s_t
# step j:  P <- (N-j+1) P ;  T <- j T + (N-j+1) P_old ;  D <- j D
# a segment (j1, j2] acts as (alpha, delta, tau):  P'=aP, T'=dT+tau P, D'=dD
class Stats:
    def __init__(self):
        self.leaf_steps = 0       # leaves visited in binary splitting (all chunks)
        self.bigmul_bits = 0      # sum of operand bits over big multiplications
        self.peak_bits = 0        # max live bignum bits in a chunk (approx.)
        self.chunks = 0
        self.items = {"main": 0, "lucas": 0, "ext": 0}
        self.item_bits = {"main": 0, "lucas": 0, "ext": 0}
        self.sweep_steps = 0
        self.t = {}


def bs(N, j1, j2, st):
    """exact (alpha, delta, tau) for steps j1+1..j2."""
    if j2 - j1 <= 24:
        a, d, t = 1, 1, 0
        for j in range(j1 + 1, j2 + 1):
            u = N - j + 1
            t = j * t + u * a
            a *= u
            d *= j
        st.leaf_steps += j2 - j1
        return mpz(a), mpz(d), mpz(t)
    mid = (j1 + j2) >> 1
    a1, d1, t1 = bs(N, j1, mid, st)
    a2, d2, t2 = bs(N, mid, j2, st)
    return a1 * a2, d1 * d2, d2 * t1 + t2 * a1


def advance(N, state, j1, j2, Q, st, lgN):
    """apply steps j1+1..j2 to state mod Q, grouping so products ~ bits(Q)."""
    P, T, D = state
    if j2 <= j1:
        return state
    qb = max(64, Q.bit_length())
    g = max(8, qb // lgN)
    j = j1
    while j < j2:
        e = min(j2, j + g)
        a, d, t = bs(N, j, e, st)
        st.bigmul_bits += 4 * (qb + a.bit_length())
        P, T, D = f_mod(a * P, Q), f_mod(d * T + t * P, Q), f_mod(d * D, Q)
        j = e
    return (P, T, D)


def art_chunk(N, items, st, lgN, out):
    """items: sorted list of (t, q, tag). Writes out[tag] = (P,T,D) mod q at t."""
    L = len(items)
    # product tree over moduli
    tree = {}

    def build(lo, hi):
        if hi - lo == 1:
            tree[(lo, hi)] = mpz(items[lo][1])
        else:
            mid = (lo + hi) >> 1
            tree[(lo, hi)] = build(lo, mid) * build(mid, hi)
        return tree[(lo, hi)]

    Q = build(0, L)
    st.peak_bits = max(st.peak_bits, Q.bit_length() * (4 + max(1, L).bit_length()))
    one = (mpz(1), mpz(1), mpz(1))
    X = advance(N, one, 0, items[0][0], Q, st, lgN)

    def rec(lo, hi, X):
        if hi - lo == 1:
            t, q, tag = items[lo]
            out[tag] = tuple(int(f_mod(x, q)) for x in X)
            return
        mid = (lo + hi) >> 1
        QL, QR = tree[(lo, mid)], tree[(mid, hi)]
        rec(lo, mid, tuple(f_mod(x, QL) for x in X))
        XR = tuple(f_mod(x, QR) for x in X)
        rec(mid, hi, advance(N, XR, items[lo][0], items[mid][0], QR, st, lgN))

    rec(0, L, X)


# ----------------------------------------------------------------- Lucas: s_t mod p, p <= t, p || m_k
def lucas_s(N, t, p, sr, cr, cache):
    """s_t mod p from s_r, C(N,r) mod p (r = t mod p) via Lucas' theorem."""
    Nd, td = [], []
    x, y = N, t
    while x:
        Nd.append(x % p)
        td.append(y % p)
        x //= p
        y //= p
    # G_i(x) = sum_{y<x} C(N_i, y) mod p ; binomial C(N_i, x) mod p
    def row(i):
        key = (p, i)
        if key not in cache:
            ni = Nd[i]
            c, pref, cs = 1, [0], []
            for yv in range(0, ni + 1):
                if yv:
                    c = c * (ni - yv + 1) % p * pow(yv, -1, p) % p
                cs.append(c)
                pref.append((pref[-1] + c) % p)
            cache[key] = (cs, pref)
        return cache[key]

    D = len(Nd)
    pow2low = [1] * (D + 1)  # prod_{u<i} 2^{N_u}
    for i in range(D):
        pow2low[i + 1] = pow2low[i] * pow(2, Nd[i], p) % p
    hi_prod = 1  # prod_{u>i} C(N_u, t_u)
    acc = 0
    for i in range(D - 1, -1, -1):
        if i == 0:
            G = (sr - cr) % p
            Cb = cr
        else:
            cs, pref = row(i)
            ti = td[i]
            G = pref[min(ti, len(cs))]
            Cb = cs[ti] if ti < len(cs) else 0
        acc = (acc + hi_prod * G % p * pow2low[i]) % p
        hi_prod = hi_prod * Cb % p
        if hi_prod == 0:
            break
    return (acc + hi_prod) % p


# ----------------------------------------------------------------- p-adic Lucas recursion for p^e, e>=2
class PadicBinom:
    """s_t(N) = sum_{j<=t} C(N,j) and C(N,t) modulo p^e for a fixed prime p.

    Uses (1+x)^N = (1+x)^{N0} ((1+x^p) + p g(x))^K,  N = pK + N0,  g = ((1+x)^p - 1 - x^p)/p,
    truncated binomial expansion  ((1+x^p)+pg)^K == sum_{i<e} C(K,i) p^i g^i (1+x^p)^{K-i}  (mod p^e),
    so coefficients/prefix sums of row N reduce to rows K-i (~N/p) modulo p^{e-i}.
    Only the low-degree polynomials h_i = (1+x)^{N0} g^i (degree < (i+1)p) are touched
    directly: cost O(e^2 p) per level and target, recursion depth log_p N.
    For e = 1 this is exactly Lucas' theorem."""

    def __init__(self, p, emax):
        self.p = p
        self.emax = emax
        self.P = [p ** i for i in range(emax + 1)]
        Pe = self.P[emax]
        # g_i = C(p,i)/p = C(p-1,i-1)/i  (i < p invertible), kept mod p^emax
        g = [0] * p
        c = 1  # C(p-1, i-1)
        for i in range(1, p):
            g[i] = c * pow(i, -1, Pe) % Pe
            c = c * (p - i) % Pe * pow(i, -1, Pe) % Pe
        gp = [[1]]  # powers of g, exact enough mod p^emax
        for i in range(1, emax):
            prev = gp[-1]
            nxt = [0] * (len(prev) + p - 1)
            for a_, x in enumerate(prev):
                if x:
                    for b_ in range(1, p):
                        nxt[a_ + b_] = (nxt[a_ + b_] + x * g[b_]) % Pe
            gp.append(nxt)
        self.gp = gp
        self.memoS, self.memoC, self.rows = {}, {}, {}
        self.work = 0

    def _row(self, N0, mod):
        """coefficients and prefix sums of (1+x)^N0 mod `mod`, N0 < p."""
        key = (N0, mod)
        if key not in self.rows:
            a, c = [], 1
            for y in range(N0 + 1):
                if y:
                    c = c * (N0 - y + 1) % mod * pow(y, -1, mod) % mod
                a.append(c)
            A, s = [], 0
            for x in a:
                s = (s + x) % mod
                A.append(s)
            self.rows[key] = (a, A)
            self.work += N0 + 1
        return self.rows[key]

    def _h(self, N0, i, u, mod, prefix):
        """(prefix sum up to u of) coefficient u of (1+x)^N0 g^i, mod `mod`."""
        a, A = self._row(N0, mod)
        gi = self.gp[i]
        tot = 0
        # h_i[u] = sum_z gi[z] a[u-z];  H_i[u] = sum_z gi[z] A[u-z]
        lo = max(0, u - N0) if not prefix else 0
        for z in range(lo, min(u, len(gi) - 1) + 1):
            if gi[z]:
                w = u - z
                if prefix:
                    tot += gi[z] * A[min(w, N0)]
                else:
                    tot += gi[z] * a[w]
        self.work += len(gi)
        return tot % mod

    def S(self, N, t, e):
        if e == 0 or t < 0:
            return 0
        mod = self.P[e]
        if t >= N:
            return pow(2, N, mod)
        key = (N, t, e)
        if key in self.memoS:
            return self.memoS[key]
        p = self.p
        if N < p:
            r = self._row(N, mod)[1][t]
        else:
            K, N0 = divmod(N, p)
            T, rr = divmod(t, p)
            r = 0
            for i in range(e):
                if K < i:
                    break
                coef = self.P[i] * math.comb(K, i)
                if coef % mod == 0:
                    continue
                m2 = self.P[e - i]
                full = pow(2, N0, m2) * (sum(self.gp[i]) % m2) % m2
                part = full * self.S(K - i, T - i - 1, e - i)
                for d in range(i + 1):
                    part += self._h(N0, i, d * p + rr, m2, True) * self.C(K - i, T - d, e - i)
                r += coef * (part % m2)
            r %= mod
        self.memoS[key] = r
        return r

    def C(self, N, t, e):
        if e == 0 or t < 0 or t > N:
            return 0
        mod = self.P[e]
        key = (N, t, e)
        if key in self.memoC:
            return self.memoC[key]
        p = self.p
        if N < p:
            r = self._row(N, mod)[0][t]
        else:
            K, N0 = divmod(N, p)
            T, rr = divmod(t, p)
            r = 0
            for i in range(e):
                if K < i:
                    break
                coef = self.P[i] * math.comb(K, i)
                if coef % mod == 0:
                    continue
                m2 = self.P[e - i]
                part = 0
                for d in range(i + 1):
                    part += self._h(N0, i, d * p + rr, m2, False) * self.C(K - i, T - d, e - i)
                r += coef * (part % m2)
            r %= mod
        self.memoC[key] = r
        return r


# ----------------------------------------------------------------- sweep for small p, e>=2
def sweep(N, p, E, targets, st):
    """s_t mod p^E for all t in targets (sorted), one pass, Gourdon Alg.2 with one prime."""
    q = p ** E
    A, B, S, v = 1, 1, 1, 0
    res = {}
    ti = 0
    tl = sorted(targets)
    if tl and tl[0] == 0:
        res[0] = 1
    for j in range(1, tl[-1] + 1):
        a, b = N - j + 1, j
        while a % p == 0:
            a //= p
            v += 1
        while b % p == 0:
            b //= p
            v -= 1
        A = A * a % q
        B = B * b % q
        S = (S * b + (A * pow(p, v, q) if v < E else 0)) % q
        if j in targets:
            res[j] = S * pow(B, -1, q) % q
    st.sweep_steps += tl[-1]
    return res


# ----------------------------------------------------------------- the C part, Theorem-2 style
def c_part_art(n, M, N, m_bits, Ps=None, st=None, fac=None, small="padic"):
    """small = 'padic': every p^e (e>=2, p<=t) part via PadicBinom (default);
       small = 'sweep': p<=Ps via O(N) sweeps, p>Ps via extended-modulus ART items."""
    st = st or Stats()
    t0 = time.perf_counter()
    if fac is None:
        fac = factor_interval(M, N)
    st.t["factor"] = time.perf_counter() - t0
    if Ps is None:
        Ps = max(3, int(math.isqrt(N)))
    base = 2 * M * N + 1
    lgN = N.bit_length()
    two_N = {}
    acc = 0
    Xk = [0] * N
    items = []          # (t, q, tag)
    lucas_need = {}     # (p, r) -> list of (k, t)
    sweep_need = {}     # p -> list of (k, t, e)
    ext_info = {}
    t0 = time.perf_counter()
    for k in range(N):
        mk = base + 2 * k
        Xk[k] = pow(5, N - 2, mk) * pow(10, n - N + 2, mk) % mk
        t = k if k <= N - 1 - k else N - 1 - k
        good = 1
        for p, e in fac[k]:
            if p > t:
                good *= p ** e
            elif e == 1:
                lucas_need.setdefault((p, t % p), []).append((k, t))
            elif small == "padic" or p <= Ps:
                sweep_need.setdefault(p, []).append((k, t, e))
            else:
                v = vp_fact(t, p)
                tag = ("ext", k, p)
                ext_info[tag] = (k, t, p, e, v)
                items.append((t, p ** (e + v), tag))
        if good > 1:
            items.append((t, good, ("main", k, good)))
    for (p, r) in lucas_need:
        items.append((r, p, ("lucas", p, r)))
    items.sort(key=lambda it: it[0])
    for it in items:
        st.items[it[2][0]] += 1
        st.item_bits[it[2][0]] += it[1].bit_length()
    st.t["setup"] = time.perf_counter() - t0

    def add(k, q, s_mod_q):
        """add (-1)^k frac(X s_k / m_k) restricted to part q (partial fractions)."""
        nonlocal acc
        mk = base + 2 * k
        if k != (k if k <= N - 1 - k else N - 1 - k):  # k was mirrored
            s_mod_q = (pow(2, N, q) - s_mod_q) % q
        a = pow((mk // q) % q, -1, q)
        y = Xk[k] % q * s_mod_q % q * a % q
        term = (y << FRAC_BITS) // q
        acc = (acc + term) % ONE if k % 2 == 0 else (acc - term) % ONE

    # --- ART over chunks
    t0 = time.perf_counter()
    out = {}
    lucas_val = {}
    i = 0
    while i < len(items):
        bits, j = 0, i
        while j < len(items) and (bits < m_bits or j == i):
            bits += items[j][1].bit_length()
            j += 1
        chunk = items[i:j]
        out.clear()
        art_chunk(N, chunk, st, lgN, out)
        st.chunks += 1
        for tag, (P, T, D) in out.items():
            kind = tag[0]
            if kind == "main":
                _, k, q = tag
                add(k, q, T * pow(D, -1, q) % q)
            elif kind == "lucas":
                _, p, r = tag
                Di = pow(D, -1, p)
                lucas_val[(p, r)] = (T * Di % p, P * Di % p)
            else:
                k, t, p, e, v = ext_info[tag]
                pe = p ** e
                Tn = (T // p ** v) % pe
                Dn = (D // p ** v) % pe
                add(k, pe, Tn * pow(Dn, -1, pe) % pe)
        i = j
    st.t["art"] = time.perf_counter() - t0

    # --- Lucas consumers
    t0 = time.perf_counter()
    cache = {}
    for (p, r), lst in lucas_need.items():
        sr, cr = lucas_val[(p, r)]
        for k, t in lst:
            add(k, p, lucas_s(N, t, p, sr, cr, cache))
    st.t["lucas"] = time.perf_counter() - t0

    # --- sweeps
    t0 = time.perf_counter()
    for p, lst in sweep_need.items():
        E = max(e for _, _, e in lst)
        if small == "padic":
            pb = PadicBinom(p, E)
            for k, t, e in lst:
                add(k, p ** e, pb.S(N, t, e))
            st.sweep_steps += pb.work
        else:
            res = sweep(N, p, E, {t for _, t, _ in lst}, st)
            for k, t, e in lst:
                add(k, p ** e, res[t] % p ** e)
    st.t["sweep"] = time.perf_counter() - t0
    return acc, st


# ----------------------------------------------------------------- Theorem-1 baseline (Algorithm 2 per term)
def c_part_thm1(n, M, N, fac=None, st=None):
    """Gourdon Algorithm 2 per k (with the k > N/2 symmetry). Returns (acc, word_op_count)."""
    if fac is None:
        fac = factor_interval(M, N)
    base = 2 * M * N + 1
    acc = 0
    ops = 0
    for k in range(N):
        m = base + 2 * k
        t = k if k <= N - 1 - k else N - 1 - k
        ps = [p for p, e in fac[k] if p <= t]
        A = B = C = 1
        R = {p: 1 for p in ps}
        Rprod = 1
        for j in range(1, t + 1):
            a, b = N - j + 1, j
            for p in ps:
                while a % p == 0:
                    a //= p
                    R[p] *= p
                    ops += 1
                while b % p == 0:
                    b //= p
                    R[p] //= p
                    ops += 1
            if ps:
                Rprod = 1
                for p in ps:
                    Rprod = Rprod * R[p] % m
            A = A * a % m
            B = B * b % m
            C = (C * b + A * Rprod) % m
            ops += 3 + 2 * len(ps)
        s = C * pow(B, -1, m) % m
        if t != k:
            s = (pow(2, N, m) - s) % m
        y = pow(5, N - 2, m) * pow(10, n - N + 2, m) % m * s % m
        term = (y << FRAC_BITS) // m
        acc = (acc + term) % ONE if k % 2 == 0 else (acc - term) % ONE
    return acc, ops


# ----------------------------------------------------------------- reference
def reference(n, digs=20):
    from mpmath import mp, mpf, floor
    mp.dps = n + digs + 20
    x = mp.pi * mpf(10) ** n
    return int((x - floor(x)) * mpf(2) ** FRAC_BITS)


def fx(v, d=15):
    return f"{v / ONE:.{d}f}"


def run(n, m_bits, M=None, n0=12, baseline=False, check=True, Ps=None, small="padic"):
    if M is None:
        # balance B (~(M+1)N powmods) against C (~N^2 lgN / m_bits leaf steps)
        M = max(4, 2 * round(n * 2.0 / m_bits))
    N = choose_N(n, n0, M)
    st = Stats()
    t0 = time.perf_counter()
    b = b_part(n, M, N)
    tB = time.perf_counter() - t0
    t0 = time.perf_counter()
    c, st = c_part_art(n, M, N, m_bits, Ps=Ps, st=st, small=small)
    tC = time.perf_counter() - t0
    val = (b - c) % ONE
    res = dict(n=n, m_bits=m_bits, M=M, N=N, tB=tB, tC=tC, val=val, st=st)
    if baseline:
        t0 = time.perf_counter()
        c1, ops = c_part_thm1(n, M, N)
        res.update(tC1=time.perf_counter() - t0, ops1=ops, agree=abs(((c1 - c + ONE // 2) % ONE) - ONE // 2) < (ONE >> 60))
    if check:
        ref = reference(n)
        res["err"] = abs(((val - ref + ONE // 2) % ONE) - ONE // 2) / ONE
    return res


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    mb = int(sys.argv[2]) if len(sys.argv) > 2 else 1024
    bl = len(sys.argv) > 3 and sys.argv[3] == "base"
    r = run(n, mb, baseline=bl)
    st = r["st"]
    print(f"n={n} m_bits={mb} M={r['M']} N={r['N']}  frac={fx(r['val'])} err={r.get('err'):.2e}")
    print(f"  tB={r['tB']:.2f}s tC={r['tC']:.2f}s  parts={ {k: round(v, 2) for k, v in st.t.items()} }")
    print(f"  chunks={st.chunks} leaf_steps={st.leaf_steps} sweep_steps={st.sweep_steps} peak_bits~{st.peak_bits}")
    print(f"  items={st.items} bits={st.item_bits}")
    if bl:
        print(f"  thm1: tC1={r['tC1']:.2f}s ops={r['ops1']} agree={r['agree']}")
