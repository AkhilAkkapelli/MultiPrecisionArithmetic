#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>
#include <stdalign.h>

size_t iamax_scalar_i32(const int32_t *restrict x, size_t n) {
    if (n == 0) return 0;
    int32_t vx = x[0];
    int32_t s = vx >> 31;
    uint32_t maxv = (uint32_t)((vx ^ s) - s);
    size_t idx = 0;

    for (size_t i = 1; i < n; i++) {
        int32_t v = x[i];
        int32_t m = v >> 31;
        uint32_t a = (uint32_t)((v ^ m) - m);
        if (a > maxv) {
            maxv = a;
            idx = i;
        }
    }
    return idx;
}

size_t iamax_avx512_i32(const int32_t *restrict x, size_t n) {
    if (n == 0) return 0;
    __m512i vmax = _mm512_setzero_si512();
    size_t i = 0;

    for (; i + 16 <= n; i += 16) {
        __m512i v = _mm512_load_si512((const __m512i *)&x[i]);
        __m512i vabs = _mm512_abs_epi32(v);
        vmax = _mm512_max_epu32(vmax, vabs);
    }

    uint32_t maxv = (uint32_t)_mm512_reduce_max_epu32(vmax);

    for (; i < n; i++) {
        int32_t v = x[i];
        int32_t s = v >> 31;
        uint32_t a = (uint32_t)((v ^ s) - s);
        if (a > maxv) maxv = a;
    }

    __m512i vmaxv = _mm512_set1_epi32((int32_t)maxv);
    for (i = 0; i + 16 <= n; i += 16) {
        __m512i v = _mm512_load_si512((const __m512i *)&x[i]);
        __m512i vabs = _mm512_abs_epi32(v);
        __mmask16 m = _mm512_cmpeq_epu32_mask(vabs, vmaxv);
        if (m) return i + (size_t)_tzcnt_u32(m);
    }

    for (; i < n; i++) {
        int32_t v = x[i];
        int32_t s = v >> 31;
        uint32_t a = (uint32_t)((v ^ s) - s);
        if (a == maxv) return i;
    }

    return 0;
}

size_t iamax_scalar_i64(const int64_t *restrict x, size_t n) {
    if (n == 0) return 0;
    int64_t vx = x[0];
    int64_t s = vx >> 63;
    uint64_t maxv = (uint64_t)((vx ^ s) - s);
    size_t idx = 0;

    for (size_t i = 1; i < n; i++) {
        int64_t v = x[i];
        int64_t m = v >> 63;
        uint64_t a = (uint64_t)((v ^ m) - m);
        if (a > maxv) {
            maxv = a;
            idx = i;
        }
    }
    return idx;
}

size_t iamax_avx512_i64(const int64_t *restrict x, size_t n) {
    if (n == 0) return 0;
    __m512i vmax = _mm512_setzero_si512();
    size_t i = 0;

    for (; i + 8 <= n; i += 8) {
        __m512i v = _mm512_load_si512((const __m512i *)&x[i]);
        __m512i vabs = _mm512_abs_epi64(v);
        vmax = _mm512_max_epu64(vmax, vabs);
    }

    uint64_t maxv = (uint64_t)_mm512_reduce_max_epu64(vmax);

    for (; i < n; i++) {
        int64_t v = x[i];
        int64_t s = v >> 63;
        uint64_t a = (uint64_t)((v ^ s) - s);
        if (a > maxv) maxv = a;
    }

    __m512i vmaxv = _mm512_set1_epi64((int64_t)maxv);
    for (i = 0; i + 8 <= n; i += 8) {
        __m512i v = _mm512_load_si512((const __m512i *)&x[i]);
        __m512i vabs = _mm512_abs_epi64(v);
        __mmask8 m = _mm512_cmpeq_epu64_mask(vabs, vmaxv);
        if (m) return i + (size_t)_tzcnt_u32(m);
    }

    for (; i < n; i++) {
        int64_t v = x[i];
        int64_t s = v >> 63;
        uint64_t a = (uint64_t)((v ^ s) - s);
        if (a == maxv) return i;
    }

    return 0;
}
