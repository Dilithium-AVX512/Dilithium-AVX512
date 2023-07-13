#ifndef CPUCYCLES_H
#define CPUCYCLES_H

#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include <stdio.h>

#ifdef USE_RDPMC /* Needs echo 2 > /sys/devices/cpu/rdpmc */

static inline uint64_t cpucycles(void)
{
  const uint32_t ecx = (1U << 30) + 1;
  uint64_t result;

  __asm__ volatile("rdpmc; shlq $32,%%rdx; orq %%rdx,%%rax"
                   : "=a"(result)
                   : "c"(ecx)
                   : "rdx");

  return result;
}

#else

static inline uint64_t cpucycles(void)
{
  uint64_t result;

  __asm__ volatile("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax"
                   : "=a"(result)
                   :
                   : "%rdx");

  return result;
}

#endif

static inline uint64_t time_now()
{
  struct timeval tv;
  uint64_t ret;

  gettimeofday(&tv, NULL);
  ret = tv.tv_sec;
  ret *= 1000000;
  ret += tv.tv_usec;

  return ret;
}

#define cycles_now cpucycles

#define oper_second_n(RANDOM, LABEL, FUNCTION, NUMBER, PARALLEL) \
  do                                                             \
  {                                                              \
    printf("%-30s : ", #LABEL);                                  \
    RANDOM;                                                      \
                                                                 \
    unsigned i;                                                  \
    uint64_t start, end;                                         \
    const unsigned iterations = NUMBER;                          \
    uint64_t start_c, end_c;                                     \
                                                                 \
    /* Load the caches*/                                         \
    for (i = 0; i < NUMBER / 10; ++i)                            \
    {                                                            \
      FUNCTION;                                                  \
    }                                                            \
                                                                 \
    start = time_now();                                          \
    start_c = cycles_now();                                      \
    for (i = 0; i < iterations; ++i)                             \
    {                                                            \
      FUNCTION;                                                  \
    }                                                            \
    end = time_now();                                            \
    end_c = cycles_now();                                        \
                                                                 \
    printf("%3lu µs, %8.1f oper/s, %6lu cycles/op\n",            \
           (unsigned long)((end - start) / iterations),          \
           PARALLEL *iterations *(double)1e6 / (end - start),    \
           (unsigned long)((end_c - start_c) / iterations));     \
  } while (0)

uint64_t cpucycles_overhead(void);

#endif