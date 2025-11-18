#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
#include <math.h>

int32_t norm_scalar_i32(const int32_t *restrict x, size_t n) {
    if (n == 0) return 0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        double v = (double)x[i];
        sum += v * v;
    }
    return (int32_t)sqrt(sum);
}

int32_t norm_avx512_i32(const int32_t *restrict x, size_t n) {
    if (n == 0) return 0;
    __m512d vsum = _mm512_setzero_pd();
    size_t i = 0;

    for (; i + 16 <= n; i += 16) {
        __m512i vi = _mm512_load_si512((const __m512i *)&x[i]);
        __m256i lo = _mm512_castsi512_si256(vi);
        __m256i hi = _mm512_extracti64x4_epi64(vi, 1);
        __m512d vlo = _mm512_cvtepi32_pd(lo);
        __m512d vhi = _mm512_cvtepi32_pd(hi);
        vsum = _mm512_fmadd_pd(vlo, vlo, vsum);
        vsum = _mm512_fmadd_pd(vhi, vhi, vsum);
    }

    double sum = _mm512_reduce_add_pd(vsum);

    for (; i < n; i++) {
        double v = (double)x[i];
        sum += v * v;
    }

    return (int32_t)sqrt(sum);
}

int64_t norm_scalar_i64(const int64_t *restrict x, size_t n) {
    if (n == 0) return 0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        double v = (double)x[i];
        sum += v * v;
    }
    return (int64_t)sqrt(sum);
}

int64_t norm_avx512_i64(const int64_t *restrict x, size_t n) {
    if (n == 0) return 0;
    __m512d vsum = _mm512_setzero_pd();
    size_t i = 0;

    for (; i + 8 <= n; i += 8) {
        __m512i vi = _mm512_load_si512((const __m512i *)&x[i]);
        __m512d vd  = _mm512_cvtepi64_pd(vi);
        vsum = _mm512_fmadd_pd(vd, vd, vsum);
    }

    double sum = _mm512_reduce_add_pd(vsum);

    for (; i < n; i++) {
        double v = (double)x[i];
        sum += v * v;
    }

    return (int64_t)sqrt(sum);
}
