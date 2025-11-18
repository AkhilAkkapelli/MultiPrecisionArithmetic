#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>

void axpy_scalar_i32(int32_t *restrict y, const int32_t *restrict x,
                     size_t n, const int32_t alpha)
{
    for (size_t i = 0; i < n; i++)
        y[i] += alpha * x[i];
}

void axpy_avx512_i32(int32_t *restrict y, const int32_t *restrict x,
                     size_t n, const int32_t alpha)
{
    __m512i alpha_vec = _mm512_set1_epi32(alpha);
    size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m512i x_vec = _mm512_load_si512((const __m512i *)&x[i]);
        __m512i y_vec = _mm512_load_si512((const __m512i *)&y[i]);
        __m512i prod  = _mm512_mullo_epi32(alpha_vec, x_vec);
        __m512i sum   = _mm512_add_epi32(y_vec, prod);
        _mm512_store_si512((__m512i *)&y[i], sum);
    }
    for (; i < n; i++)
        y[i] += alpha * x[i];
}

void axpy_scalar_i64(int64_t *restrict y, const int64_t *restrict x,
                     size_t n, const int64_t alpha)
{
    for (size_t i = 0; i < n; i++)
        y[i] += alpha * x[i];
}

void axpy_avx512_i64(int64_t *restrict y, const int64_t *restrict x,
                     size_t n, const int64_t alpha)
{
    __m512i alpha_vec = _mm512_set1_epi64(alpha);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i x_vec = _mm512_load_si512((const __m512i *)&x[i]);
        __m512i y_vec = _mm512_load_si512((const __m512i *)&y[i]);
        __m512i prod  = _mm512_mullo_epi64(alpha_vec, x_vec);
        __m512i sum   = _mm512_add_epi64(y_vec, prod);
        _mm512_store_si512((__m512i *)&y[i], sum);
    }
    for (; i < n; i++)
        y[i] += alpha * x[i];
}
