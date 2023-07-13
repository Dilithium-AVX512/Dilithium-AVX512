#include <stdint.h>
#include "../sign.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../params.h"
#include "cpucycles.h"
#include "speed_print.h"
#include "../rejsample.h"


#define NTESTS 1000000

uint64_t t[NTESTS];

int main(void)
{
//   unsigned int i;
  size_t siglen;
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t sig[CRYPTO_BYTES];
  uint8_t seed[CRHBYTES];
  polyvecl mat[K];
  poly *a = &mat[0].vec[0];
  poly *b = &mat[0].vec[1];
  poly *c = &mat[0].vec[2];
  

  uint8_t aaa[848];
  uint8_t *bbb=aaa;

  oper_second_n(while (0), rej_uniform_avx, rej_uniform_avx(a->coeffs, bbb), NTESTS, 1);
  oper_second_n(while (0), rej_eta_avx, rej_eta_avx(a->coeffs, bbb), NTESTS, 1);
  oper_second_n(while (0), polyvec_matrix_expand, polyvec_matrix_expand(mat, seed), NTESTS, 1);
  oper_second_n(while (0), poly_uniform_eta, poly_uniform_eta(a, seed, 0), NTESTS, 1);
  oper_second_n(while (0), poly_uniform_gamma1, poly_uniform_gamma1(a, seed, 0), NTESTS, 1);
  oper_second_n(while (0), poly_ntt, poly_ntt(a), NTESTS, 1);
  oper_second_n(while (0), poly_invntt_tomont, poly_invntt_tomont(a), NTESTS, 1);
  oper_second_n(while (0), poly_pointwise_montgomery, poly_pointwise_montgomery(c, a, b), NTESTS, 1);
  oper_second_n(while (0), poly_challenge, poly_challenge(c, seed), NTESTS, 1);
  oper_second_n(while (0), crypto_sign_keypair, crypto_sign_keypair(pk, sk), NTESTS, 1);
  oper_second_n(while (0), crypto_sign_signature, crypto_sign_signature(sig, &siglen, sig, CRHBYTES, sk), NTESTS, 1);
  oper_second_n(while (0), crypto_sign_verify, crypto_sign_verify(sig, CRYPTO_BYTES, sig, CRHBYTES, pk), NTESTS, 1);


  return 0;
}
