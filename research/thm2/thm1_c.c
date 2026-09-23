/* Theorem-1 baseline for the C part only: Gourdon's Algorithm 2 per term,
 * with the k > N/2 symmetry, 64-bit word arithmetic (moduli < 2^32).
 * plus the B part, i.e. the complete Theorem-1 algorithm (Gourdon's M = 2*ceil(n/log^3 n)).
 * usage: thm1_c n M N    -> prints frac(10^n pi), frac(C) and seconds
 * build: gcc -O3 -march=native -o thm1_c thm1_c.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

typedef uint64_t u64;
typedef int64_t i64;

static u64 powmod(u64 b, u64 e, u64 m) {
    u64 r = 1 % m; b %= m;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    return r;
}
static u64 invmod(u64 a, u64 m) {
    i64 t = 0, nt = 1, r = (i64)m, nr = (i64)(a % m);
    while (nr) { i64 q = r / nr, x; x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
    if (t < 0) t += (i64)m;
    return (u64)t;
}

int main(int argc, char **argv) {
    u64 n = strtoull(argv[1], 0, 10), M = strtoull(argv[2], 0, 10), N = strtoull(argv[3], 0, 10);
    u64 base = 2 * M * N + 1;
    long double acc = 0;
    u64 ops = 0;
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (u64 k = 0; k < N; k++) {
        u64 m = base + 2 * k;
        u64 t = (k <= N - 1 - k) ? k : N - 1 - k;
        /* primes p <= t dividing m */
        u64 ps[32]; int np = 0; u64 x = m;
        for (u64 p = 3; p * p <= x; p += 2)
            if (x % p == 0) { if (p <= t) ps[np++] = p; while (x % p == 0) x /= p; }
        if (x > 1 && x <= t) ps[np++] = x;
        int v[32] = {0};
        u64 na[32], nb[32];
        for (int i = 0; i < np; i++) {
            u64 p = ps[i];
            nb[i] = p;                       /* next j with p | j */
            u64 r = (N + 1) % p;             /* j == N+1 mod p  <=> p | N-j+1 */
            na[i] = r ? r : p;
        }
        u64 A = 1, B = 1, C = 1, R = 1;
        for (u64 j = 1; j <= t; j++) {
            u64 a = N - j + 1, b = j;
            int changed = 0;
            for (int i = 0; i < np; i++) {
                u64 p = ps[i];
                if (na[i] == j) { do { a /= p; v[i]++; } while (a % p == 0); na[i] += p; changed = 1; }
                if (nb[i] == j) { do { b /= p; v[i]--; } while (b % p == 0); nb[i] += p; changed = 1; }
            }
            if (changed) { R = 1; for (int i = 0; i < np; i++) R = R * powmod(ps[i], v[i], m) % m; }
            A = A * (a % m) % m;
            B = B * (b % m) % m;
            C = (C * (b % m) + A * R) % m;
            ops++;
        }
        u64 s = C * invmod(B, m) % m;
        if (t != k) s = (powmod(2, N, m) + m - s) % m;
        u64 X = powmod(5, N - 2, m) * powmod(10, n - N + 2, m) % m;
        u64 y = X * s % m;
        long double term = (long double)y / (long double)m;
        acc += (k & 1) ? -term : term;
        acc -= (long double)(i64)acc;
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (acc < 0) acc += 1;
    /* B part: sum_{k<(M+1)N} (-1)^k (4*10^n mod 2k+1)/(2k+1) */
    struct timespec t2;
    long double bacc = 0;
    u64 K = (M + 1) * N;
    for (u64 k = 0; k < K; k++) {
        u64 q = 2 * k + 1;
        u64 x = 4 * powmod(10, n, q) % q;
        long double term = (long double)x / (long double)q;
        bacc += (k & 1) ? -term : term;
        if ((k & 1023) == 0) bacc -= (long double)(i64)bacc;
    }
    bacc -= (long double)(i64)bacc;
    clock_gettime(CLOCK_MONOTONIC, &t2);
    long double f = bacc - acc; f -= (long double)(i64)f; if (f < 0) f += 1;
    printf("frac=%.15Lf fracC=%.15Lf tC=%.3f tB=%.3f inner_steps=%llu\n", f, acc,
           (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec),
           (t2.tv_sec - t1.tv_sec) + 1e-9 * (t2.tv_nsec - t1.tv_nsec), (unsigned long long)ops);
    return 0;
}
