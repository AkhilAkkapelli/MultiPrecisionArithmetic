#define _POSIX_C_SOURCE 200809L
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

extern void axpy_scalar(int64_t *y, const int64_t *x, size_t n, int64_t alpha);
extern void axpy_avx512(int64_t *y, const int64_t *x, size_t n, int64_t alpha);

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
    int64_t *x, *y1, *y2;
    int64_t alpha = 3;

    printf("\n Benchmarking C AXPY (y += alpha*x, int64) — scalar vs AVX-512\n");
    printf(" ============================================================================\n");
    printf("        N        |  Scalar(s)  |  AVX512(s)  |  Speedup  |   Check\n");
    printf(" ----------------------------------------------------------------------\n");

    for (int t = 0; t < ntests; t++) {
        size_t N = Nvals[t];
        x  = aligned_alloc(64, N * sizeof(int64_t));
        y1 = aligned_alloc(64, N * sizeof(int64_t));
        y2 = aligned_alloc(64, N * sizeof(int64_t));

        srand(0);
        for (size_t i = 0; i < N; i++) {
            x[i]  = (rand() % 2000000L) - 1000000L;
            y1[i] = (rand() % 2000000L) - 1000000L;
            y2[i] = y1[i];
        }

        double tsum_s = 0.0, tsum_v = 0.0;

        for (int rep = 0; rep < nrep; rep++) {
            memcpy(y2, y1, N * sizeof(int64_t));

            double t0 = now();
            axpy_scalar(y1, x, N, alpha);
            double t1 = now();
            tsum_s += (t1 - t0);

            t0 = now();
            axpy_avx512(y2, x, N, alpha);
            t1 = now();
            tsum_v += (t1 - t0);
        }

        tsum_s /= nrep;
        tsum_v /= nrep;

        int ok = 1;
        for (size_t i = 0; i < N; i++) {
            if (y1[i] != y2[i]) { ok = 0; break; }
        }

        printf("%12zu | %10.6f | %10.6f | %9.2fx | %s\n",
               N, tsum_s, tsum_v, tsum_s / tsum_v, ok ? "OK" : "FAIL");

        free(x);
        free(y1);
        free(y2);
    }

    printf(" ----------------------------------------------------------------------\n");
    printf(" Times averaged over %d runs.\n", nrep);
    return 0;
}
