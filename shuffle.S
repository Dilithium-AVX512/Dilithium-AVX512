#include "consts.h"
.include "shuffle.inc"

.text
nttunpack512_avx:
#load
vmovdqa32	  0+0(%rdi),%zmm4
vmovdqa32	 64+0(%rdi),%zmm5
vmovdqa32	128+0(%rdi),%zmm6
vmovdqa32	192+0(%rdi),%zmm7
vmovdqa32	256+0(%rdi),%zmm8
vmovdqa32	320+0(%rdi),%zmm9
vmovdqa32	384+0(%rdi),%zmm10
vmovdqa32	448+0(%rdi),%zmm11
vmovdqa32	512+0(%rdi),%zmm12
vmovdqa32	576+0(%rdi),%zmm13
vmovdqa32	640+0(%rdi),%zmm14
vmovdqa32	704+0(%rdi),%zmm15
vmovdqa32	768+0(%rdi),%zmm16
vmovdqa32	832+0(%rdi),%zmm17
vmovdqa32	896+0(%rdi),%zmm18
vmovdqa32	960+0(%rdi),%zmm19

shuffle16   4,12,3,12
shuffle16   5,13,4,13
shuffle16   6,14,5,14
shuffle16   7,15,6,15
shuffle16   8,16,7,16
shuffle16   9,17,8,17
shuffle16   10,18,9,18
shuffle16   11,19,10,19

shuffle8    3,7,11,7
shuffle8    12,16,3,16
shuffle8    4,8,12,8
shuffle8    13,17,4,17
shuffle8    5,9,13,9
shuffle8    14,18,5,18
shuffle8    6,10,14,10
shuffle8    15,19,6,19

shuffle4    11,13,15,13
shuffle4    7,9,11,9
shuffle4    3,5,7,5
shuffle4    16,18,3,18
shuffle4    12,14,16,14
shuffle4    8,10,12,10
shuffle4    4,6,8,6
shuffle4    17,19,4,19

shuffle2    15,16,17,16
shuffle2    13,14,15,14
shuffle2    11,12,13,12
shuffle2    9,10,11,10
shuffle2    7,8,9,8
shuffle2    5,6,7,6
shuffle2    3,4,5,4
shuffle2    18,19,3,19

#store
vmovdqa32	%zmm17,0+0(%rdi)
vmovdqa32	%zmm16,64+0(%rdi)
vmovdqa32	%zmm15,128+0(%rdi)
vmovdqa32	%zmm14,192+0(%rdi)
vmovdqa32	%zmm13,256+0(%rdi)
vmovdqa32	%zmm12,320+0(%rdi)
vmovdqa32	%zmm11,384+0(%rdi)
vmovdqa32	%zmm10,448+0(%rdi)
vmovdqa32	%zmm9,512+0(%rdi)
vmovdqa32	%zmm8,576+0(%rdi)
vmovdqa32	%zmm7,640+0(%rdi)
vmovdqa32	%zmm6,704+0(%rdi)
vmovdqa32	%zmm5,768+0(%rdi)
vmovdqa32	%zmm4,832+0(%rdi)
vmovdqa32	%zmm3,896+0(%rdi)
vmovdqa32	%zmm19,960+0(%rdi)
ret

.global cdecl(nttunpack_avx)
cdecl(nttunpack_avx):

kmovw       (_ZETAS_QINV+5)*4(%rsi),%k1 
kmovw       (_ZETAS_QINV+6)*4(%rsi),%k2
call		nttunpack512_avx

ret

