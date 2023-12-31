#ifndef NTT_H
#define NTT_H

#include <immintrin.h>

#define ntt_avx DILITHIUM_NAMESPACE(ntt_avx)
void ntt_avx(__m512i *a, const __m512i *qdata);
#define invntt_avx DILITHIUM_NAMESPACE(invntt_avx)
void invntt_avx(__m512i *a, const __m512i *qdata);

#define nttunpack_avx DILITHIUM_NAMESPACE(nttunpack_avx)
void nttunpack_avx(__m512i *a,const __m512i *qdata);

#define pointwise_avx DILITHIUM_NAMESPACE(pointwise_avx)
void pointwise_avx(__m512i *c, const __m512i *a, const __m512i *b, const __m512i *qdata);
#define pointwise_acc_avx DILITHIUM_NAMESPACE(pointwise_acc_avx)
void pointwise_acc_avx(__m512i *c, const __m512i *a, const __m512i *b, const __m512i *qdata);

#endif
