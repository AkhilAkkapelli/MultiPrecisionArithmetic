#define _POSIX_C_SOURCE 199309L
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// Correct function signatures matching the implementation
extern void scal_scalar_i64(int64_t *x, int64_t alpha, size_t n);
extern void scal_avx512_i64(int64_t *x, int64_t alpha, size_t n);

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
    const int64_t alpha = 3;

    printf("\n Benchmarking C SCAL (scale int64 vector) — scalar vs AVX-512\n");
    printf(" ============================================================================\n");
    printf("        N        |  Scalar(s)  |  AVX512(s)  |  Speedup  |   Check\n");
    printf(" ----------------------------------------------------------------------------\n");

    for (int t = 0; t < ntests; t++) {
        size_t N = Nvals[t];
        
        // Allocate aligned memory
        int64_t *a_orig   = aligned_alloc(64, N * sizeof(int64_t));
        int64_t *a_scalar = aligned_alloc(64, N * sizeof(int64_t));
        int64_t *a_avx512 = aligned_alloc(64, N * sizeof(int64_t));

        if (!a_orig || !a_scalar || !a_avx512) {
            fprintf(stderr, "Memory allocation failed\n");
            return 1;
        }

        // Initialize with random values
        srand(0);
        for (size_t i = 0; i < N; i++) {
            a_orig[i] = (rand() % 2000000L) - 1000000L;
        }

        double tsum_scalar = 0.0, tsum_avx512 = 0.0;

        // Run multiple repetitions for accurate timing
        for (int rep = 0; rep < nrep; rep++) {
            // Test scalar version
            memcpy(a_scalar, a_orig, N * sizeof(int64_t));
            double t0 = now();
            scal_scalar_i64(a_scalar, alpha, N);
            double t1 = now();
            tsum_scalar += (t1 - t0);

            // Test AVX-512 version
            memcpy(a_avx512, a_orig, N * sizeof(int64_t));
            t0 = now();
            scal_avx512_i64(a_avx512, alpha, N);
            t1 = now();
            tsum_avx512 += (t1 - t0);
        }

        // Calculate average times
        tsum_scalar /= nrep;
        tsum_avx512 /= nrep;

        // Verify results match
        int ok = 1;
        for (size_t i = 0; i < N; i++) {
            if (a_scalar[i] != a_avx512[i]) {
                ok = 0;
                fprintf(stderr, "Mismatch at index %zu: scalar=%ld, avx512=%ld\n",
                        i, a_scalar[i], a_avx512[i]);
                break;
            }
        }

        double speedup = tsum_scalar / tsum_avx512;
        printf("%12zu | %10.6f | %10.6f | %8.2fx | %s\n",
               N, tsum_scalar, tsum_avx512, speedup, ok ? "OK" : "FAIL");

        free(a_orig);
        free(a_scalar);
        free(a_avx512);
    }

    printf(" ----------------------------------------------------------------------------\n");
    printf(" Times averaged over %d runs.\n", nrep);
    printf(" Note: Speedup = Scalar_time / AVX512_time (higher is better)\n");
    
    return 0;
}
