/*
 * n-th binary (hex) digit of pi with Bellard's 1997 formula (the "pi_bin" article on
 * bellard.org/pi/pi_bin/pi_bin.html). Not Bellard's code (the page ships none): written for
 * this benchmark from the formula, in the style of the BBP digit-extraction algorithm.
 *
 *   pi = 1/2^6 * sum_{n>=0} (-1)^n / 2^(10n) * ( -2^5/(4n+1) - 1/(4n+3) + 2^8/(10n+1)
 *                                                - 2^6/(10n+3) - 2^2/(10n+5) - 2^2/(10n+7)
 *                                                + 1/(10n+9) )
 *
 * usage: pi_bin d      -> hex digits d+1 .. d+8 of pi's fractional part (d=0: 243F6A88)
 * build: gcc -O3 -march=native -o pi_bin pi_bin.c -lm
 * Time is O(d log d): it is a different (binary, much faster) problem than the decimal
 * extractors, so it is benchmarked separately and not comparable digit-for-digit with pihunt.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

static u64 pow2_mod(u64 e, u64 m) {
    u64 r = 1, b = 2 % m;
    while (e) {
        if (e & 1) r = (u128)r * b % m;
        b = (u128)b * b % m;
        e >>= 1;
    }
    return r;
}

/* frac(2^e / m) for any integer e (e may be negative). */
static long double frac2(long e, u64 m) {
    if (e >= 0) return (long double)pow2_mod((u64)e, m) / (long double)m;
    return ldexpl(1.0L, (int)e) / (long double)m;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: pi_bin d\n"); return 1; }
    long d = atol(argv[1]);
    long e = 4 * d;                       /* we want frac(2^e * pi) */
    static const struct { int j; int sign; int a, b; } T[7] = {
        { 5, -1, 4, 1}, {0, -1, 4, 3}, {8, +1, 10, 1}, {6, -1, 10, 3},
        { 2, -1, 10, 5}, {2, -1, 10, 7}, {0, +1, 10, 9}};
    long double sum = 0;
    for (long n = 0; 10 * n <= e + 40; n++) {
        long double alt = (n & 1) ? -1.0L : 1.0L, part = 0;
        for (int i = 0; i < 7; i++) {
            u64 m = (u64)T[i].a * n + T[i].b;
            part += T[i].sign * frac2(e - 6 - 10 * n + T[i].j, m);
        }
        sum += alt * part;
        sum -= floorl(sum);
    }
    sum -= floorl(sum);
    printf("Hex digits of pi at position %ld: ", d + 1);
    for (int i = 0; i < 8; i++) {
        sum *= 16;
        int h = (int)sum;
        printf("%X", h);
        sum -= h;
    }
    printf("\n");
    return 0;
}
