#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

extern uint64_t asum_scalar_i64(const int64_t *a, size_t n);
extern uint64_t asum_avx512_i64(const int64_t *a, size_t n);

static double now_seconds(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

int main(void) {
    const size_t Nvals[] = { (1u<<22), (1u<<24), (1u<<26), (1u<<28), (1u<<29) }; // 1M .. 268M
    const int ncases = (int)(sizeof(Nvals)/sizeof(Nvals[0]));
    const int trials = 5;

    printf("\n Benchmarking C ASUM (64-bit integers) — scalar vs AVX-512\n");
    printf(" ==========================================================================\n");
    printf("        N        |  Scalar(s)  |  AVX512(s)  |  AVX512x  |   Check\n");
    printf(" --------------------------------------------------------------------------\n");

    srand(12345);

    for (int c = 0; c < ncases; ++c) {
        size_t N = Nvals[c];

        // allocate aligned memory (64-byte) — check for allocation failure
        int64_t *a = (int64_t*)aligned_alloc(64, N * sizeof(int64_t));
        if (!a) {
            fprintf(stderr, "Allocation failed for N=%zu\n", N);
            return 1;
        }

        // fill with pseudo-random signed 32-bit-range values to keep sums reasonable
        for (size_t i = 0; i < N; ++i) a[i] = (int64_t)((int32_t)(rand() - RAND_MAX/2));

        double t_scalar = 0.0, t_avx = 0.0;
        uint64_t s_scalar = 0, s_avx = 0;

        for (int r = 0; r < trials; ++r) {
            double t0 = now_seconds();
            s_scalar = asum_scalar_i64(a, N);
            t_scalar += now_seconds() - t0;

            t0 = now_seconds();
            s_avx = asum_avx512_i64(a, N);
            t_avx += now_seconds() - t0;
        }
        t_scalar /= trials;
        t_avx /= trials;
        double speedup = (t_avx > 0.0) ? (t_scalar / t_avx) : 0.0;
        const char *ok = (s_scalar == s_avx) ? "OK" : "FAIL";

        printf("%12zu | %11.6f | %11.6f | %8.2fx | %7s\n",
               N, t_scalar, t_avx, speedup, ok);

        free(a);
    }

    printf(" --------------------------------------------------------------------------\n");
    printf(" Times averaged over %d runs per size.\n\n", trials);
    return 0;
}
