#define _POSIX_C_SOURCE 200809L
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

extern size_t iamin_scalar(const int64_t *x, size_t n);
extern size_t iamin_avx512(const int64_t *x, size_t n);

static inline double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

int main(void)
{
    const size_t Nvals[] = {
        1UL << 12,  // 4K
        1UL << 14,  // 16K
        1UL << 16,  // 64K
        1UL << 18,  // 256K
        1UL << 20,  // 1M
        1UL << 22,  // 4M
        1UL << 24,  // 16M
        1UL << 26,   // 64M
        1UL << 27
    };
    const int ntests = sizeof(Nvals) / sizeof(Nvals[0]);
    const int nrep = 3;  // show 3 repetitions
    int64_t *x = NULL;

    printf("\n Benchmarking C IAMIN (index of min abs value, int64)\n");
    printf(" ======================================================================================================================\n");
    printf("        N        |   Scalar(s) [3 runs]       |   AVX512(s) [3 runs]       |  Speedup[3 runs]  |  Avg(Scalar)  Avg(AVX)  AvgSpeedup  | Check\n");
    printf(" ----------------------------------------------------------------------------------------------------------------------\n");

    for (int t = 0; t < ntests; t++) {
        size_t N = Nvals[t];
        if (posix_memalign((void**)&x, 64, N * sizeof(int64_t)) != 0) {
            fprintf(stderr, "posix_memalign failed for N=%zu\n", N);
            return 1;
        }

        srand(0);
        for (size_t i = 0; i < N; i++)
            x[i] = (rand() % 2000000L) - 1000000L;

        double ts_s[3], ts_v[3], spd[3];
        size_t idx_s = 0, idx_v = 0;

        for (int rep = 0; rep < nrep; rep++) {
            double t0 = now();
            idx_s = iamin_scalar(x, N);
            double t1 = now();
            ts_s[rep] = t1 - t0;

            t0 = now();
            idx_v = iamin_avx512(x, N);
            double t2 = now();
            ts_v[rep] = t2 - t0;

            spd[rep] = (ts_v[rep] > 0.0) ? ts_s[rep] / ts_v[rep] : 0.0;
        }

        double avg_s = (ts_s[0] + ts_s[1] + ts_s[2]) / 3.0;
        double avg_v = (ts_v[0] + ts_v[1] + ts_v[2]) / 3.0;
        double avg_spd = avg_s / avg_v;
        int ok = (idx_s == idx_v);

        printf("%12zu | %8.6f %8.6f %8.6f | %8.6f %8.6f %8.6f | %6.2fx %6.2fx %6.2fx | %10.6f %10.6f %8.2fx | %s\n",
               N,
               ts_s[0], ts_s[1], ts_s[2],
               ts_v[0], ts_v[1], ts_v[2],
               spd[0], spd[1], spd[2],
               avg_s, avg_v, avg_spd,
               ok ? "OK" : "FAIL");

        free(x);
        x = NULL;
    }

    printf(" ----------------------------------------------------------------------------------------------------------------------\n");
    printf(" Shows per-run times, per-run speedups, and averages across 3 runs.\n");
    return 0;
}
