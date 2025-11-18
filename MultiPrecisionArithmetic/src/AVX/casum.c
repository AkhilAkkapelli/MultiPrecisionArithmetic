#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
#include <stdalign.h>
int32_t asum_scalar_i32(const int32_t *restrict a, size_t n) {
    int64_t sum = 0;
    for (size_t i = 0; i < n; ++i) {
        int32_t x = a[i];
        int32_t s = x >> 31;
        int32_t absx = (x ^ s) - s;
        sum += absx;
    }
    return (int32_t)sum;
}
int32_t asum_avx512_i32(const int32_t *restrict a, size_t n) {
    size_t i = 0;
    __m512i vsum = _mm512_setzero_si512();
    for (; i + 16 <= n; i += 16) {
        __m512i v = _mm512_load_si512((const __m512i *)(a + i));
        __m512i vabs = _mm512_abs_epi32(v);
        vsum = _mm512_add_epi32(vsum, vabs);
    }
    int32_t sum = _mm512_reduce_add_epi32(vsum);
    if (i < n)
        sum += asum_scalar_i32(a + i, n - i);
    return sum;
}
int64_t asum_scalar_i64(const int64_t *restrict a, size_t n) {
    int64_t sum = 0;
    for (size_t i = 0; i < n; ++i) {
        int64_t x = a[i];
        int64_t s = x >> 63;
        int64_t absx = (x ^ s) - s;
        sum += absx;
    }
    return sum;
}
int64_t asum_avx512_i64(const int64_t *restrict a, size_t n) {
    size_t i = 0;
    __m512i vsum = _mm512_setzero_si512();
    for (; i + 8 <= n; i += 8) {
        __m512i v = _mm512_load_si512((const __m512i *)(a + i));
        __m512i vabs = _mm512_abs_epi64(v);
        vsum = _mm512_add_epi64(vsum, vabs);
    }
    int64_t sum = _mm512_reduce_add_epi64(vsum);
    if (i < n)
        sum += asum_scalar_i64(a + i, n - i);
    return sum;
}
