#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>

void scal_scalar_i32(int32_t *restrict x, int32_t alpha, size_t n)
{
    if (n == 0 || alpha == 1)
        return;

    for (size_t i = 0; i < n; i++)
        x[i] = x[i] * alpha;
}

void scal_avx512_i32(int32_t *restrict x, int32_t alpha, size_t n)
{
    if (n == 0 || alpha == 1)
        return;

    const __m512i valpha = _mm512_set1_epi32(alpha);
    size_t i = 0;

    for (; i + 16 <= n; i += 16) {
        __m512i v = _mm512_load_si512((const __m512i *)&x[i]);
        v = _mm512_mullo_epi32(v, valpha);
        _mm512_store_si512((__m512i *)&x[i], v);
    }

    for (; i < n; i++)
        x[i] = x[i] * alpha;
}

void scal_scalar_i64(int64_t *restrict x, int64_t alpha, size_t n)
{
    if (n == 0 || alpha == 1)
        return;

    for (size_t i = 0; i < n; i++)
        x[i] = x[i] * alpha;
}

void scal_avx512_i64(int64_t *restrict x, int64_t alpha, size_t n)
{
    if (n == 0 || alpha == 1)
        return;

    const __m512i valpha = _mm512_set1_epi64(alpha);
    size_t i = 0;

    for (; i + 8 <= n; i += 8) {
        __m512i v = _mm512_load_si512((const __m512i *)&x[i]);
        v = _mm512_mullo_epi64(v, valpha);
        _mm512_store_si512((__m512i *)&x[i], v);
    }

    for (; i < n; i++)
        x[i] = x[i] * alpha;
}
