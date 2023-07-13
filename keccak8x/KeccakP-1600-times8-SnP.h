#ifndef _KeccakP_1600_times8_SnP_h_
#define _KeccakP_1600_times8_SnP_h_

#include "KeccakP-SIMD512-config.h"
#include "../fips202x8.h"

#define KeccakP1600times8_implementation        "512-bit SIMD implementation (" KeccakP1600times8_implementation_config ")"
#define KeccakP1600times8_statesSizeInBytes     1600
#define KeccakP1600times8_statesAlignment       64
#define KeccakF1600times8_FastLoop_supported
#define KeccakP1600times8_12rounds_FastLoop_supported
#define KeccakF1600times8_FastKravatte_supported
#define KeccakP1600times8_K12ProcessLeaves_supported

#include <stddef.h>
#include <stdint.h>

#define KeccakP1600times8_StaticInitialize()
#define KeccakP1600times8_InitializeAll FIPS202X8_NAMESPACE(KeccakP1600times8_InitializeAll)
void KeccakP1600times8_InitializeAll(void *states);
#define KeccakP1600times8_AddByte(states, instanceIndex, byte, offset) \
    ((unsigned char*)(states))[(instanceIndex)*8 + ((offset)/8)*8*8 + (offset)%8] ^= (byte)
#define KeccakP1600times8_AddBytes FIPS202X8_NAMESPACE(KeccakP1600times8_AddBytes)
void KeccakP1600times8_AddBytes(void *states, unsigned int instanceIndex, const unsigned char *data, unsigned int offset, unsigned int length);
#define KeccakP1600times8_AddLanesAll FIPS202X8_NAMESPACE(KeccakP1600times8_AddLanesAll)
void KeccakP1600times8_AddLanesAll(void *states, const unsigned char *data, unsigned int laneCount, unsigned int laneOffset);
#define KeccakP1600times8_OverwriteBytes FIPS202X8_NAMESPACE(KeccakP1600times8_OverwriteBytes)
void KeccakP1600times8_OverwriteBytes(void *states, unsigned int instanceIndex, const unsigned char *data, unsigned int offset, unsigned int length);
#define KeccakP1600times8_OverwriteLanesAll FIPS202X8_NAMESPACE(KeccakP1600times8_OverwriteLanesAll)
void KeccakP1600times8_OverwriteLanesAll(void *states, const unsigned char *data, unsigned int laneCount, unsigned int laneOffset);
#define KeccakP1600times8_OverwriteWithZeroes FIPS202X8_NAMESPACE(KeccakP1600times8_OverwriteWithZeroes)
void KeccakP1600times8_OverwriteWithZeroes(void *states, unsigned int instanceIndex, unsigned int byteCount);
#define KeccakP1600times8_PermuteAll_4rounds FIPS202X8_NAMESPACE(KeccakP1600times8_PermuteAll_4rounds)
void KeccakP1600times8_PermuteAll_4rounds(void *states);
#define KeccakP1600times8_PermuteAll_6rounds FIPS202X8_NAMESPACE(KeccakP1600times8_PermuteAll_6rounds)
void KeccakP1600times8_PermuteAll_6rounds(void *states);
#define KeccakP1600times8_PermuteAll_12rounds FIPS202X8_NAMESPACE(KeccakP1600times8_PermuteAll_12rounds)
void KeccakP1600times8_PermuteAll_12rounds(void *states);
#define KeccakP1600times8_PermuteAll_24rounds FIPS202X8_NAMESPACE(KeccakP1600times8_PermuteAll_24rounds)
void KeccakP1600times8_PermuteAll_24rounds(void *states);
#define KeccakP1600times8_ExtractBytes FIPS202X8_NAMESPACE(KeccakP1600times8_ExtractBytes)
void KeccakP1600times8_ExtractBytes(const void *states, unsigned int instanceIndex, unsigned char *data, unsigned int offset, unsigned int length);
#define KeccakP1600times8_ExtractLanesAll FIPS202X8_NAMESPACE(KeccakP1600times8_ExtractLanesAll)
void KeccakP1600times8_ExtractLanesAll(const void *states, unsigned char *data, unsigned int laneCount, unsigned int laneOffset);
#define KeccakP1600times8_ExtractAndAddBytes FIPS202X8_NAMESPACE(KeccakP1600times8_ExtractAndAddBytes)
void KeccakP1600times8_ExtractAndAddBytes(const void *states, unsigned int instanceIndex,  const unsigned char *input, unsigned char *output, unsigned int offset, unsigned int length);
#define KeccakP1600times8_ExtractAndAddLanesAll FIPS202X8_NAMESPACE(KeccakP1600times8_ExtractAndAddLanesAll)
void KeccakP1600times8_ExtractAndAddLanesAll(const void *states, const unsigned char *input, unsigned char *output, unsigned int laneCount, unsigned int laneOffset);
#define KeccakF1600times8_FastLoop_Absorb FIPS202X8_NAMESPACE(KeccakF1600times8_FastLoop_Absorb)
size_t KeccakF1600times8_FastLoop_Absorb(void *states, unsigned int laneCount, unsigned int laneOffsetParallel, unsigned int laneOffsetSerial, const unsigned char *data, size_t dataByteLen);
#define KeccakP1600times8_12rounds_FastLoop_Absorb FIPS202X8_NAMESPACE(KeccakP1600times8_12rounds_FastLoop_Absorb)
size_t KeccakP1600times8_12rounds_FastLoop_Absorb(void *states, unsigned int laneCount, unsigned int laneOffsetParallel, unsigned int laneOffsetSerial, const unsigned char *data, size_t dataByteLen);
#define KeccakP1600times8_KravatteCompress FIPS202X8_NAMESPACE(KeccakP1600times8_KravatteCompress)
size_t KeccakP1600times8_KravatteCompress(uint64_t *xAccu, uint64_t *kRoll, const unsigned char *input, size_t inputByteLen);
#define KeccakP1600times8_KravatteExpand FIPS202X8_NAMESPACE(KeccakP1600times8_KravatteExpand)
size_t KeccakP1600times8_KravatteExpand(uint64_t *yAccu, const uint64_t *kRoll, unsigned char *output, size_t outputByteLen);
#define KeccakP1600times8_K12ProcessLeaves FIPS202X8_NAMESPACE(KeccakP1600times8_K12ProcessLeaves)
void KeccakP1600times8_K12ProcessLeaves(const unsigned char *input, unsigned char *output);

#endif
