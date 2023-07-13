#include <stdint.h>
#include <immintrin.h>
#include <string.h>
#include "params.h"
#include "rounding.h"
#include "rejsample.h"
#include "consts.h"

#define _mm512_blendv_epi32(a, b, zero, f) \
  _mm512_mask_blend_epi32(_mm512_cmpgt_epi32_mask(zero, f), a, b)

/*************************************************
* Name:        power2round
*
* Description: For finite field elements a, compute a0, a1 such that
*              a mod^+ Q = a1*2^D + a0 with -2^{D-1} < a0 <= 2^{D-1}.
*              Assumes a to be positive standard representative.
*
* Arguments:   - __m256i *a1: output array of length N/8 with high bits
*              - __m256i *a0: output array of length N/8 with low bits a0
*              - const __m256i *a: input array of length N/8
*
**************************************************/
void power2round_avx(__m512i *a1, __m512i *a0, const __m512i *a)
{
  unsigned int i;
  __m512i f, f0, f1;
  const __m512i mask = _mm512_set1_epi32(-(1 << D));
  const __m512i half = _mm512_set1_epi32((1 << (D - 1)) - 1);

  for (i = 0; i < N / 16; ++i)
  {
    f = _mm512_load_si512(&a[i]);
    f1 = _mm512_add_epi32(f, half);
    f0 = _mm512_and_si512(f1, mask);
    f1 = _mm512_srli_epi32(f1, D);
    f0 = _mm512_sub_epi32(f, f0);
    _mm512_store_si512(&a1[i], f1);
    _mm512_store_si512(&a0[i], f0);
  }
}

/*************************************************
* Name:        decompose
*
* Description: For finite field element a, compute high and low parts a0, a1 such
*              that a mod^+ Q = a1*ALPHA + a0 with -ALPHA/2 < a0 <= ALPHA/2 except
*              if a1 = (Q-1)/ALPHA where we set a1 = 0 and
*              -ALPHA/2 <= a0 = a mod Q - Q < 0. Assumes a to be positive standard
*              representative.
*
* Arguments:   - __m256i *a1: output array of length N/8 with high parts
*              - __m256i *a0: output array of length N/8 with low parts a0
*              - const __m256i *a: input array of length N/8
*
**************************************************/
#if GAMMA2 == (Q - 1) / 32
void decompose_avx(__m512i *a1, __m512i *a0, const __m512i *a)
{
  unsigned int i;
  __m512i f, f0, f1;
  __mmask16 mmask;
  const __m512i q = _mm512_load_epi32(&qdata.vec512[_16XQ / 16]);
  const __m512i hq = _mm512_srli_epi32(q, 1);
  const __m512i v = _mm512_set1_epi32(1025);
  const __m512i alpha = _mm512_set1_epi32(2 * GAMMA2);
  const __m512i off = _mm512_set1_epi32(127);
  const __m512i shift = _mm512_set1_epi32(512);
  const __m512i mask = _mm512_set1_epi32(15);
  const __m512i zero=_mm512_setzero_si512();

  for (i = 0; i < N / 16; i++)
  {
    f = _mm512_load_si512(&a[i]);
    f1 = _mm512_add_epi32(f, off);
    f1 = _mm512_srli_epi32(f1, 7);
    f1 = _mm512_mulhi_epu16(f1, v);
    f1 = _mm512_mulhrs_epi16(f1, shift);
    f1 = _mm512_and_si512(f1, mask);
    f0 = _mm512_mullo_epi32(f1, alpha);
    f0 = _mm512_sub_epi32(f, f0);
    mmask = _mm512_cmpgt_epi32_mask(f0, hq);
    f = _mm512_mask_blend_epi32(mmask, zero, q);
    f0 = _mm512_sub_epi32(f0, f);
    _mm512_store_si512(&a1[i], f1);
    _mm512_store_si512(&a0[i], f0);
  }
}

#elif GAMMA2 == (Q - 1) / 88
void decompose_avx(__m512i *a1, __m512i *a0, const __m512i *a)
{
  unsigned int i;
  __m512i f, f0, f1, t;
  __mmask16 mmask;
  const __m512i q = _mm512_load_si512(&qdata.vec512[_16XQ / 16]);
  const __m512i hq = _mm512_srli_epi32(q, 1);
  const __m512i v = _mm512_set1_epi32(11275);
  const __m512i alpha = _mm512_set1_epi32(2 * GAMMA2);
  const __m512i off = _mm512_set1_epi32(127);
  const __m512i shift = _mm512_set1_epi32(128);
  const __m512i max = _mm512_set1_epi32(43);
  const __m512i zero = _mm512_setzero_si512();

  for (i = 0; i < N / 16; i++)
  {
    f = _mm512_load_si512(&a[i]);
    f1 = _mm512_add_epi32(f, off);
    f1 = _mm512_srli_epi32(f1, 7);
    f1 = _mm512_mulhi_epu16(f1, v);
    f1 = _mm512_mulhrs_epi16(f1, shift);
    t = _mm512_sub_epi32(max, f1);
    f1 = _mm512_blendv_epi32(f1, zero, zero, t);
    f0 = _mm512_mullo_epi32(f1, alpha);
    f0 = _mm512_sub_epi32(f, f0);
    mmask = _mm512_cmpgt_epi32_mask(f0, hq);
    f = _mm512_mask_blend_epi32(mmask, _mm512_setzero_si512(), q);
    f0 = _mm512_sub_epi32(f0, f);
    _mm512_store_si512(&a1[i], f1);
    _mm512_store_si512(&a0[i], f0);
  }
}
#endif

/*************************************************
* Name:        make_hint
*
* Description: Compute indices of polynomial coefficients whose low bits
*              overflow into the high bits.
*
* Arguments:   - uint8_t *hint: hint array
*              - const __m256i *a0: low bits of input elements
*              - const __m256i *a1: high bits of input elements
*
* Returns number of overflowing low bits
**************************************************/
unsigned int make_hint_avx(uint8_t hint[N], const __m512i *restrict a0, const __m512i *restrict a1)
{
  unsigned int i, n = 0;
  __m512i f0, f1, g0, g1;
  uint32_t bad;
  uint64_t idx;
  __m256i t;
  __mmask16 mmask;
  const __m512i low = _mm512_set1_epi32(-GAMMA2);
  const __m512i high = _mm512_set1_epi32(GAMMA2);

  for (i = 0; i < N / 16; ++i)
  {
    f0 = _mm512_load_si512(&a0[i]);
    f1 = _mm512_load_si512(&a1[i]);
    g0 = _mm512_abs_epi32(f0);
    mmask = _mm512_cmpgt_epi32_mask(g0, high);
    g0 = _mm512_mask_set1_epi32(_mm512_setzero_si512(), mmask, 0XFFFFFFFF);
    mmask = _mm512_cmpeq_epi32_mask(f0, low);
    g1 = _mm512_mask_blend_epi32(mmask, _mm512_set1_epi32(0), _mm512_xor_epi32(f1,_mm512_set1_epi32(-1)));
    g1 = _mm512_mask_blend_epi32(_mm512_cmpeq_epi32_mask(f1, _mm512_setzero_si512()),g1, _mm512_setzero_si512());
    g0 = _mm512_or_si512(g0, g1);

    t = _mm512_castsi512_si256(g0);
    bad = _mm256_movemask_ps((__m256)t);
    memcpy(&idx, idxlut[bad], 8);
    idx += (uint64_t)0x0808080808080808 * (2*i);
    memcpy(&hint[n], &idx, 8);
    n += _mm_popcnt_u32(bad);

    t = _mm512_extracti32x8_epi32(g0, 1);
    bad = _mm256_movemask_ps((__m256)t);
    memcpy(&idx, idxlut[bad], 8);
    idx += (uint64_t)0x0808080808080808 * (2*i + 1);
    memcpy(&hint[n], &idx, 8);
    n += _mm_popcnt_u32(bad);
  }

  return n;
}


/*************************************************
* Name:        use_hint
*
* Description: Correct high parts according to hint.
*
* Arguments:   - __m256i *b: output array of length N/8 with corrected high parts
*              - const __m256i *a: input array of length N/8
*              - const __m256i *a: input array of length N/8 with hint bits
*
**************************************************/
void use_hint_avx(__m512i *b, const __m512i *a, const __m512i *restrict hint)
{
  unsigned int i;
  __m512i a0[N / 16];
  __m512i f, g, h, t;
  
  const __m512i zero = _mm512_setzero_si512();
#if GAMMA2 == (Q - 1) / 32
  const __m512i mask = _mm512_set1_epi32(15);
#elif GAMMA2 == (Q - 1) / 88
  __mmask16 mmask;
  const __m512i max = _mm512_set1_epi32(43);
#endif

  decompose_avx(b, a0, a);
  for (i = 0; i < N / 16; i++)
  {
    f = _mm512_load_si512(&a0[i]);
    g = _mm512_load_si512(&b[i]);
    h = _mm512_load_si512(&hint[i]);
    t = _mm512_blendv_epi32(zero, h, zero, f);
    t = _mm512_slli_epi32(t, 1);
    h = _mm512_sub_epi32(h, t);
    g = _mm512_add_epi32(g, h);
#if GAMMA2 == (Q - 1) / 32
    g = _mm512_and_si512(g, mask);
#elif GAMMA2 == (Q - 1) / 88
    g = _mm512_blendv_epi32(g, max, zero, g);
    mmask = _mm512_cmpgt_epi32_mask(g, max);
    g = _mm512_mask_blend_epi32(mmask, g, zero);
#endif
    _mm512_store_si512(&b[i], g);
  }
}

