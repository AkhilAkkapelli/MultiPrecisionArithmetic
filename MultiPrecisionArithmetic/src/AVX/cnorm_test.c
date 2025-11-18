#define _POSIX_C_SOURCE 200809L
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

extern int64_t norm_scalar_i64(const int64_t *a, size_t n);
extern int64_t norm_avx512_i64(const int64_t *a, size_t n);

static double now()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

int main()
{
    const size_t Nvals[] = { 1UL<<20, 1UL<<22, 1UL<<24, 1UL<<26, 1UL<<28 };
    const int ntests = sizeof(Nvals)/sizeof(Nvals[0]);
    const int nrep = 5;

    int64_t *a;
    printf("\n Benchmarking C NORM (sum of squares, int64) — scalar vs AVX-512\n");
    printf(" ============================================================================\n");
    printf("        N        |  Scalar(s)  |  AVX512(s)  |  Speedup  |   Check\n");
    printf(" ----------------------------------------------------------------------\n");

    for (int t = 0; t < ntests; t++) {
        size_t N = Nvals[t];
        a = aligned_alloc(64, N * sizeof(int64_t));

        srand(0);
        for (size_t i = 0; i < N; i++)
            a[i] = (rand() % 2000000L) - 1000000L;

        double tsum_s = 0.0, tsum_v = 0.0;
        int64_t ref = 0, res = 0;

        for (int rep = 0; rep < nrep; rep++) {
            double t0 = now();
            ref = norm_scalar_i64(a, N);
            double t1 = now();
            tsum_s += (t1 - t0);

            t0 = now();
            res = norm_avx512_i64(a, N);
            t1 = now();
            tsum_v += (t1 - t0);
        }

        tsum_s /= nrep;
        tsum_v /= nrep;

        printf("%12zu | %10.6f | %10.6f | %9.2fx | %s\n",
               N, tsum_s, tsum_v, tsum_s / tsum_v, (ref == res) ? "OK" : "FAIL");

        free(a);
    }

    printf(" ----------------------------------------------------------------------\n");
    printf(" Times averaged over %d runs.\n", nrep);
    return 0;
}
