#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <immintrin.h>
#include <string.h>
#include "fips202.h"
#include "fips202x8.h"

/* Use implementation from the Keccak Code Package */
#define KeccakF1600_StatePermute8x FIPS202X8_NAMESPACE(KeccakP1600times8_PermuteAll_24rounds)
extern void KeccakF1600_StatePermute8x(__m512i *s);

static void keccakx8_absorb_once(__m512i s[25],
                                 unsigned int r,
                                 const uint8_t *in0,
                                 const uint8_t *in1,
                                 const uint8_t *in2,
                                 const uint8_t *in3,
                                 const uint8_t *in4,
                                 const uint8_t *in5,
                                 const uint8_t *in6,
                                 const uint8_t *in7,
                                 size_t inlen,
                                 uint8_t p)
{
  size_t i;
  uint64_t pos = 0;
  __m512i t, idx;

  for(i = 0; i < 25; ++i)
    s[i] = _mm512_setzero_si512();

/*
0001
0002
__m512i _mm512_i64gather_epi64 (__m512i vindex, void const* base_addr=8, int scale=1)
FOR j := 0 to 7
	i := j*64
	m := j*64
	addr := base_addr + vindex[m+63:m] * ZeroExtend64(scale) * 8
  addr = 8+in0*8
	dst[i+63:i] := MEM[addr+63:addr]
  dst[63:0]=MEM[in0*8+63:in0*8]
ENDFOR
dst[MAX:512] := 0
*/
  idx = _mm512_set_epi64((long long)in7, (long long)in6, (long long)in5, (long long)in4,(long long)in3, (long long)in2, (long long)in1, (long long)in0);
  while(inlen >= r) {
    for(i = 0; i < r/8; ++i) {
      t = _mm512_i64gather_epi64(idx, (long long *)pos, 1);
      s[i] = _mm512_xor_si512(s[i], t);
      pos += 8;
    }
    inlen -= r;

    KeccakF1600_StatePermute8x(s);
  }

  // uint64_t *pp, *pp2;
  // uint8_t *bb;
  // pp=in0;
  // pp2=in0+8;
  // bb=in0;
  // printf("in0, in0+8 %lu %lu\n", (*pp), (*pp2));
  for(i = 0; i < inlen/8; ++i) {
    t = _mm512_i64gather_epi64(idx, (long long *)pos, 1);
    // pp = (uint64_t*)&t;
    // printf("t %lu\n", *pp);
    s[i] = _mm512_xor_si512(s[i], t);
    pos += 8;
  }
  inlen -= 8*i;

  if(inlen) {
    t = _mm512_i64gather_epi64(idx, (long long *)pos, 1);
    idx = _mm512_set1_epi64((1ULL << (8*inlen)) - 1);
    t = _mm512_and_si512(t, idx);
    s[i] = _mm512_xor_si512(s[i], t);
  }

  t = _mm512_set1_epi64((uint64_t)p << 8*inlen);
  s[i] = _mm512_xor_si512(s[i], t);
  t = _mm512_set1_epi64(1ULL << 63);
  s[r/8 - 1] = _mm512_xor_si512(s[r/8 - 1], t);
}

static void keccakx8_squeezeblocks(uint8_t *out0,
                                   uint8_t *out1,
                                   uint8_t *out2,
                                   uint8_t *out3,
                                   uint8_t *out4,
                                   uint8_t *out5,
                                   uint8_t *out6,
                                   uint8_t *out7,
                                   size_t nblocks,
                                   unsigned int r,
                                   __m512i s[25])
{
  unsigned int i;
  __m128d t;

  while(nblocks > 0) {
    KeccakF1600_StatePermute8x(s);
    for(i=0; i < r/8; ++i) {
      t = _mm_castsi128_pd(_mm512_castsi512_si128(s[i]));
      _mm_storel_pd((__attribute__((__may_alias__)) double *)&out0[8*i], t);
      _mm_storeh_pd((__attribute__((__may_alias__)) double *)&out1[8*i], t);
      t = _mm_castsi128_pd(_mm512_extracti64x2_epi64 (s[i],1));
      _mm_storel_pd((__attribute__((__may_alias__)) double *)&out2[8*i], t);
      _mm_storeh_pd((__attribute__((__may_alias__)) double *)&out3[8*i], t);
      t = _mm_castsi128_pd(_mm512_extracti64x2_epi64 (s[i],2));
      _mm_storel_pd((__attribute__((__may_alias__)) double *)&out4[8*i], t);
      _mm_storeh_pd((__attribute__((__may_alias__)) double *)&out5[8*i], t);
      t = _mm_castsi128_pd(_mm512_extracti64x2_epi64 (s[i],3));
      _mm_storel_pd((__attribute__((__may_alias__)) double *)&out6[8*i], t);
      _mm_storeh_pd((__attribute__((__may_alias__)) double *)&out7[8*i], t);
    }

    out0 += r;
    out1 += r;
    out2 += r;
    out3 += r;
    out4 += r;
    out5 += r;
    out6 += r;
    out7 += r;
    --nblocks;
  }
}

void shake128x8_absorb_once(keccakx8_state *state,
                            const uint8_t *in0,
                            const uint8_t *in1,
                            const uint8_t *in2,
                            const uint8_t *in3,
                            const uint8_t *in4,
                            const uint8_t *in5,
                            const uint8_t *in6,
                            const uint8_t *in7,
                            size_t inlen)
{
  keccakx8_absorb_once(state->s, SHAKE128_RATE, in0, in1, in2, in3, in4, in5, in6, in7, inlen, 0x1F);
}

void shake128x8_squeezeblocks(uint8_t *out0,
                              uint8_t *out1,
                              uint8_t *out2,
                              uint8_t *out3,
                              uint8_t *out4,
                              uint8_t *out5,
                              uint8_t *out6,
                              uint8_t *out7,
                              size_t nblocks,
                              keccakx8_state *state)
{
  keccakx8_squeezeblocks(out0, out1, out2, out3, out4, out5, out6, out7, nblocks, SHAKE128_RATE, state->s);
}

void shake256x8_absorb_once(keccakx8_state *state,
                            const uint8_t *in0,
                            const uint8_t *in1,
                            const uint8_t *in2,
                            const uint8_t *in3,
                            const uint8_t *in4,
                            const uint8_t *in5,
                            const uint8_t *in6,
                            const uint8_t *in7,
                            size_t inlen)
{
  keccakx8_absorb_once(state->s, SHAKE256_RATE, in0, in1, in2, in3, in4, in5, in6, in7, inlen, 0x1F);
}

void shake256x8_squeezeblocks(uint8_t *out0,
                              uint8_t *out1,
                              uint8_t *out2,
                              uint8_t *out3,
                              uint8_t *out4,
                              uint8_t *out5,
                              uint8_t *out6,
                              uint8_t *out7,
                              size_t nblocks,
                              keccakx8_state *state)
{
  keccakx8_squeezeblocks(out0, out1, out2, out3, out4, out5, out6, out7, nblocks, SHAKE256_RATE, state->s);
}

void shake128x8(uint8_t *out0,
                uint8_t *out1,
                uint8_t *out2,
                uint8_t *out3,
                uint8_t *out4,
                uint8_t *out5,
                uint8_t *out6,
                uint8_t *out7,
                size_t outlen,
                const uint8_t *in0,
                const uint8_t *in1,
                const uint8_t *in2,
                const uint8_t *in3,
                const uint8_t *in4,
                const uint8_t *in5,
                const uint8_t *in6,
                const uint8_t *in7,
                size_t inlen)
{
  unsigned int i;
  size_t nblocks = outlen/SHAKE128_RATE;
  uint8_t t[8][SHAKE128_RATE];
  keccakx8_state state;

  shake128x8_absorb_once(&state, in0, in1, in2, in3, in4, in5, in6, in7, inlen);
  shake128x8_squeezeblocks(out0, out1, out2, out3, out4, out5, out6, out7, nblocks, &state);

  out0 += nblocks*SHAKE128_RATE;
  out1 += nblocks*SHAKE128_RATE;
  out2 += nblocks*SHAKE128_RATE;
  out3 += nblocks*SHAKE128_RATE;
  out4 += nblocks*SHAKE128_RATE;
  out5 += nblocks*SHAKE128_RATE;
  out6 += nblocks*SHAKE128_RATE;
  out7 += nblocks*SHAKE128_RATE;
  outlen -= nblocks*SHAKE128_RATE;

  if(outlen) {
    shake128x8_squeezeblocks(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], 1, &state);
    for(i = 0; i < outlen; ++i) {
      out0[i] = t[0][i];
      out1[i] = t[1][i];
      out2[i] = t[2][i];
      out3[i] = t[3][i];
      out4[i] = t[4][i];
      out5[i] = t[5][i];
      out6[i] = t[6][i];
      out7[i] = t[7][i];
    }
  }
}

void shake256x8(uint8_t *out0,
                uint8_t *out1,
                uint8_t *out2,
                uint8_t *out3,
                uint8_t *out4,
                uint8_t *out5,
                uint8_t *out6,
                uint8_t *out7,
                size_t outlen,
                const uint8_t *in0,
                const uint8_t *in1,
                const uint8_t *in2,
                const uint8_t *in3,
                const uint8_t *in4,
                const uint8_t *in5,
                const uint8_t *in6,
                const uint8_t *in7,
                size_t inlen)
{
  unsigned int i;
  size_t nblocks = outlen/SHAKE256_RATE;
  uint8_t t[8][SHAKE256_RATE];
  keccakx8_state state;

  shake256x8_absorb_once(&state, in0, in1, in2, in3, in4, in5, in6, in7, inlen);
  shake256x8_squeezeblocks(out0, out1, out2, out3, out4, out5, out6, out7, nblocks, &state);

  out0 += nblocks*SHAKE256_RATE;
  out1 += nblocks*SHAKE256_RATE;
  out2 += nblocks*SHAKE256_RATE;
  out3 += nblocks*SHAKE256_RATE;
  out4 += nblocks*SHAKE256_RATE;
  out5 += nblocks*SHAKE256_RATE;
  out6 += nblocks*SHAKE256_RATE;
  out7 += nblocks*SHAKE256_RATE;
  outlen -= nblocks*SHAKE256_RATE;

  if(outlen) {
    shake256x8_squeezeblocks(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], 1, &state);
    for(i = 0; i < outlen; ++i) {
      out0[i] = t[0][i];
      out1[i] = t[1][i];
      out2[i] = t[2][i];
      out3[i] = t[3][i];
      out4[i] = t[4][i];
      out5[i] = t[5][i];
      out6[i] = t[6][i];
      out7[i] = t[7][i];
    }
  }
}


// int main(){
//   keccakx8_state state;
//   uint8_t buf[8][848] = {71,213,108,32,174,108,248,27,246,61,185,102,13,170,206,220,111,188,218,246,144,71,168,133,212,163,203,174,230,33,0,106,1,248,21,181,252,16,175,25,123,32,90,135,35,57,29,250,76,201,194,190,203,39,57,131,174,190,53,190,28,219,125,22,109,141,175,172,118,133,88,158,146,215,84,252,171,54,125,45,249,45,145,245,99,125,85,222,220,240,23,22,71,216,46,222,38,246,94,252,74,221,227,216,153,151,25,144,1,63,230,141,237,142,169,8,228,251,3,125,55,73,200,65,225,202,12,42,88,171,106,108,2,86,2,135,103,51,48,165,115,87,136,158,33,114,230,86,2,206,74,230,59,60,103,29,67,134,83,215,61,9,180,176,150,202,232,6,88,211,159,182,196,47,135,99,53,72,159,202,213,41,17,51,153,80,90,55,170,122,45,126,178,214,90,45,67,140,163,241,51,161,52,178,33,14,254,59,214,15,125,149,22,198,160,201,123,131,151,95,230,31,248,57,114,131,36,214,181,131,235,129,252,25,200,157,22,166,21,160,62,148,246,212,212,3,4,202,71,240,210,220,33,40,12,147,99,38,180,172,225,228,57,41,201,58,78,182,75,8,218,65,6,158,220,133,192,82,13,80,99,150,177,26,48,252,23,176,163,195,114,143,221,30,68,166,154,192,62,15,212,235,185,23,76,47,101,34,252,168,61,246,130,199,41,73,204,27,50,46,226,172,178,195,185,25,5,20,214,173,57,152,120,89,167,6,53,50,49,74,116,124,142,78,154,34,111,59,223,138,17,122,167,76,224,193,201,172,245,203,119,251,201,254,139,167,67,246,15,154,76,74,29,172,234,123,214,135,233,139,255,15,0,56,206,144,56,65,219,249,37,200,228,121,138,209,5,188,53,155,49,16,247,4,64,127,84,154,175,142,86,115,142,123,186,45,163,124,155,74,225,252,150,123,160,255,134,179,143,123,40,124,182,86,57,16,17,115,176,146,102,173,68,125,213,21,80,229,240,241,33,2,160,145,184,33,114,39,205,240,165,46,222,193,166,194,76,19,170,165,159,243,184,125,6,107,243,79,213,196,122,234,141,116,151,236,124,229,97,171,181,12,233,81,144,181,254,121,22,232,125,58,174,88,28,103,116,110,88,210,3,245,200,196,37,49,79,215,155,119,254,161,112,240,31,14,80,247,212,139,220,185,211,75,69,188,217,146,118,26,105,144,235,182,75,106,98,244,253,65,1,201,155,104,36,202,121,177,120,241,14,180,128,231,87,159,241,248,224,192,252,172,29,156,220,114,19,130,163,150,210,150,113,131,117,219,216,156,144,193,60,109,114,39,135,42,46,219,133,75,55,244,50,75,91,67,85,14,195,76,175,70,190,249,212,123,177,106,245,55,76,77,133,148,93,20,249,192,183,134,158,205,19,119,131,191,17,94,24,238,169,120,211,52,106,76,189,122,56,12,112,54,154,190,123,187,158,160,57,125,153,118,178,79,158,253,186,247,7,225,239,85,200,141,79,128,46,146,152,7,254,31,157,14,211,168,106,8,150,57,222,185,214,23,6,47,181,147,227,148,88,123,154,150,219,76,243,199,90,107,126,51,99,80,103,199,44,153,115,116,239,11,211,196,200,173,81,25,103,82,59,182,240,107,143,62,153,44,244,80,166,213,105,207,85,103,22,175,81,96,148,250,185,226,33,149,22,179,133,165,243,177,78,166,242,239,197,61,58,2,126,90,87,138,125,194,194,2,108,228,175,228,207,147,205,241,184,124,104,134,202,18,20,59,109,3,204,189,21,82,244,86,132,195,128,230,128,149,134,67,24,186,216,252,252,105,10,142,198,122,174,120,61,200,131,149,189,107,35,227,70,114,235,250,170,239,40,196,148,203,0,0,0,0,0,0,0,0,
//                          73,223,217,128,155,188,84,1,74,171,204,106,154,25,245,237,72,173,87,217,25,2,145,114,1,182,137,120,42,198,199,94,109,169,170,16,46,52,46,199,95,77,25,100,85,184,200,185,188,117,95,153,247,89,22,106,234,151,26,138,188,192,199,21,239,201,61,158,82,139,105,70,86,224,17,107,96,55,117,160,213,149,170,178,65,149,250,105,47,101,130,161,197,198,134,44,244,9,210,197,2,177,133,220,245,31,95,130,85,91,36,164,57,137,91,217,130,71,71,189,23,198,198,243,71,146,13,249,247,58,95,74,8,85,80,153,159,137,124,81,155,58,205,146,124,107,145,212,226,75,234,189,29,133,13,197,240,222,72,6,15,194,21,32,126,67,38,225,72,3,165,118,191,31,169,16,76,234,19,231,93,251,7,9,211,115,244,143,19,172,56,240,31,53,100,186,73,87,220,92,45,161,193,132,99,115,89,1,28,254,157,98,83,200,112,232,162,205,82,35,250,114,174,182,223,173,34,200,5,33,87,175,241,184,165,143,90,18,97,155,118,1,105,44,116,0,122,67,187,217,71,89,188,16,20,92,18,247,191,7,151,104,65,253,149,226,76,7,11,24,61,225,148,170,75,165,173,223,47,177,119,176,58,174,80,3,73,36,75,129,217,21,10,187,159,27,150,19,92,136,160,253,189,28,183,203,152,99,102,189,46,113,26,146,108,199,75,131,41,169,106,44,10,172,128,52,31,221,179,255,41,226,213,84,172,19,52,159,151,212,249,126,199,80,127,151,173,182,51,210,236,142,250,9,115,217,176,27,89,197,199,83,215,163,136,149,253,96,120,126,149,94,121,29,172,236,101,134,188,208,249,29,78,166,54,63,84,8,248,24,110,90,117,58,25,72,242,228,254,145,238,208,18,240,83,74,59,176,38,200,84,61,31,16,104,41,63,100,125,168,136,33,142,37,230,140,44,237,157,77,175,205,247,13,4,138,205,96,147,139,7,123,12,181,46,189,83,161,217,102,94,108,153,199,55,123,29,119,52,101,157,107,11,46,112,223,163,122,189,119,201,193,58,22,19,152,220,233,59,182,58,116,190,56,41,75,60,199,192,204,14,227,206,146,248,7,245,28,69,145,146,65,168,82,17,77,2,96,176,78,89,157,228,14,106,183,9,249,133,73,88,56,67,7,173,238,77,123,15,54,175,210,140,51,212,251,199,195,236,86,35,56,69,209,139,254,18,193,162,0,212,251,32,171,240,213,199,58,92,168,130,123,21,9,143,104,63,238,64,46,172,198,77,102,135,14,238,244,223,135,115,214,120,8,239,238,6,180,140,168,237,62,159,218,86,50,31,239,191,244,84,125,186,242,24,151,162,81,24,1,190,19,88,190,245,69,74,108,230,56,102,12,50,220,233,195,76,253,164,103,249,70,51,62,47,172,91,74,198,92,201,145,84,84,184,208,192,193,216,248,10,189,229,255,21,89,31,243,63,24,240,71,217,164,191,58,221,56,119,231,31,238,165,252,112,167,214,60,118,143,252,247,33,56,17,58,7,146,206,158,129,48,11,116,49,21,234,210,100,87,39,109,128,121,44,224,185,59,134,73,68,111,197,175,116,245,32,147,202,178,204,147,44,207,208,81,209,32,30,217,143,26,197,203,91,234,140,183,109,19,192,137,1,177,200,7,211,208,241,108,66,63,95,248,121,190,55,119,65,65,115,146,33,68,163,217,120,144,228,142,131,8,83,242,27,83,150,76,4,206,66,40,200,207,17,117,80,246,6,188,188,37,253,184,153,125,205,193,138,81,7,88,44,152,15,99,0,130,30,250,7,223,51,241,39,80,48,77,230,191,4,146,181,182,39,136,115,192,15,50,36,152,42,25,64,34,71,54,64,55,0,0,0,0,0,0,0,0};

//    for(int i=0;i<8;i++){
//     for(int j=0;j<848;j++){
//       printf("%u ",buf[i][j]);
//     }
//     printf("\n");
//   }


//   shake128x8_absorb_once(&state,buf[0],buf[1],buf[2],buf[3],buf[4],buf[5],buf[6],buf[7],34);
//   shake128x8_squeezeblocks(buf[0],buf[1],buf[2],buf[3],buf[4],buf[5],buf[6],buf[7],5,&state);
//   for(int i=0;i<8;i++){
//     for(int j=0;j<848;j++){
//       printf("%u ",buf[i][j]);
//     }
//     printf("\n");
//   }
//   return 0;

// }


