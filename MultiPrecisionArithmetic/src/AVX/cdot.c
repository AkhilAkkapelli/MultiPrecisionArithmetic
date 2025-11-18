#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
#include <stdalign.h>

int32_t dot_scalar_i32(const int32_t *restrict a, const int32_t *restrict b, size_t n)
{
    int64_t sum = 0;
    for (size_t i = 0; i < n; i++)
        sum += a[i] * b[i];
    return (int32_t)sum;
}

int32_t dot_avx512_i32(const int32_t *restrict a, const int32_t *restrict b, size_t n)
{
    __m512i vsum = _mm512_setzero_si512();
    size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m512i va = _mm512_load_si512((const __m512i *)&a[i]);
        __m512i vb = _mm512_load_si512((const __m512i *)&b[i]);
        __m512i vmul = _mm512_mullo_epi32(va, vb);
        vsum = _mm512_add_epi32(vsum, vmul);
    }
    int32_t sum = _mm512_reduce_add_epi32(vsum);
    for (; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}

int64_t dot_scalar_i64(const int64_t *restrict a, const int64_t *restrict b, size_t n)
{
    int64_t sum = 0;
    for (size_t i = 0; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}

int64_t dot_avx512_i64(const int64_t *restrict a, const int64_t *restrict b, size_t n)
{
    __m512i vsum = _mm512_setzero_si512();
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i va = _mm512_load_si512((const __m512i *)&a[i]);
        __m512i vb = _mm512_load_si512((const __m512i *)&b[i]);
        __m512i vmul = _mm512_mullo_epi64(va, vb);
        vsum = _mm512_add_epi64(vsum, vmul);
    }
    int64_t sum = _mm512_reduce_add_epi64(vsum);
    for (; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}
