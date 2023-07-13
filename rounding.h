#ifndef ROUNDING_H
#define ROUNDING_H

#include <stdint.h>
#include <immintrin.h>
#include "params.h"

#define power2round_avx DILITHIUM_NAMESPACE(power2round_avx)
void power2round_avx(__m512i *a1, __m512i *a0, const __m512i *a);
#define decompose_avx DILITHIUM_NAMESPACE(decompose_avx)
void decompose_avx(__m512i *a1, __m512i *a0, const __m512i *a);
#define make_hint_avx DILITHIUM_NAMESPACE(make_hint_avx)
unsigned int make_hint_avx(uint8_t hint[N], const __m512i *a0, const __m512i *a1);
#define use_hint_avx DILITHIUM_NAMESPACE(use_hint_avx)
void use_hint_avx(__m512i *b, const __m512i *a, const __m512i *hint);

#endif
