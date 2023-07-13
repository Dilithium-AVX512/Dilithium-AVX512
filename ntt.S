#include "consts.h"
.include "shuffle.inc"

.macro butterfly l,h,zl0=1,zl1=1,zh0=2,zh1=2
vpmuldq		%zmm\zl0,%zmm\h,%zmm29
vmovshdup	%zmm\h,%zmm28
vpmuldq		%zmm\zl1,%zmm28,%zmm30

vpmuldq		%zmm\zh0,%zmm\h,%zmm\h
vpmuldq		%zmm\zh1,%zmm28,%zmm28

vpmuldq		%zmm0,%zmm29,%zmm29
vpmuldq		%zmm0,%zmm30,%zmm30

vmovshdup	%zmm\h,%zmm\h 
vpblendmd	%zmm28,%zmm\h,%zmm\h {%k1}

vpsubd		%zmm\h,%zmm\l,%zmm28
vpaddd		%zmm\h,%zmm\l,%zmm\l

vmovshdup	%zmm29,%zmm29
vpblendmd	%zmm30,%zmm29,%zmm29 {%k1}

vpaddd		%zmm29,%zmm28,%zmm\h
vpsubd		%zmm29,%zmm\l,%zmm\l
.endm

.macro levels0t7
/* level 0 */
vmovdqa32		  0+0(%rdi),%zmm4
vmovdqa32		 64+0(%rdi),%zmm5
vmovdqa32		128+0(%rdi),%zmm6
vmovdqa32	 	192+0(%rdi),%zmm7
vmovdqa32		256+0(%rdi),%zmm8
vmovdqa32		320+0(%rdi),%zmm9
vmovdqa32		384+0(%rdi),%zmm10
vmovdqa32	 	448+0(%rdi),%zmm11
vmovdqa32		512+0(%rdi),%zmm12
vmovdqa32		576+0(%rdi),%zmm13
vmovdqa32		640+0(%rdi),%zmm14
vmovdqa32	 	704+0(%rdi),%zmm15
vmovdqa32		768+0(%rdi),%zmm16
vmovdqa32		832+0(%rdi),%zmm17
vmovdqa32		896+0(%rdi),%zmm18
vmovdqa32	 	960+0(%rdi),%zmm19

vpbroadcastd	(_ZETAS_QINV+7)*4(%rsi),%zmm1
vpbroadcastd	(_ZETAS+7)*4(%rsi),%zmm2
butterfly	4,12
butterfly	5,13
butterfly	6,14
butterfly	7,15
butterfly	8,16
butterfly	9,17
butterfly	10,18
butterfly	11,19

shuffle16   4,12,3,12
shuffle16   5,13,4,13
shuffle16   6,14,5,14
shuffle16   7,15,6,15
shuffle16   8,16,7,16
shuffle16   9,17,8,17
shuffle16   10,18,9,18
shuffle16   11,19,10,19


/*level 1*/
vmovdqa32     (_ZETAS_QINV+16)*4(%rsi),%zmm1 
vmovdqa32     (_ZETAS+16)*4(%rsi),%zmm2   

butterfly   3,7
butterfly	12,16
butterfly	4,8
butterfly	13,17
butterfly	5,9
butterfly	14,18
butterfly	6,10
butterfly	15,19

shuffle8    3,7,11,7
shuffle8    12,16,3,16
shuffle8    4,8,12,8
shuffle8    13,17,4,17
shuffle8    5,9,13,9
shuffle8    14,18,5,18
shuffle8    6,10,14,10
shuffle8    15,19,6,19


/* level 2 */
vmovdqa32 	(_ZETAS_QINV+32)*4(%rsi),%zmm1
vmovdqa32 	(_ZETAS+32)*4(%rsi),%zmm2

butterfly	11,13
butterfly	7,9
butterfly	3,5
butterfly	16,18
butterfly	12,14
butterfly	8,10
butterfly	4,6
butterfly	17,19

shuffle4    11,13,15,13
shuffle4    7,9,11,9
shuffle4    3,5,7,5
shuffle4    16,18,3,18
shuffle4    12,14,16,14
shuffle4    8,10,12,10
shuffle4    4,6,8,6
shuffle4    17,19,4,19


/* level 3 */
vmovdqa32	(_ZETAS_QINV+48)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+48)*4(%rsi),%zmm2

butterfly	15,16
butterfly	13,14
butterfly	11,12
butterfly	9,10
butterfly	7,8
butterfly	5,6
butterfly	3,4
butterfly	18,19

shuffle2    15,16,17,16
shuffle2    13,14,15,14
shuffle2    11,12,13,12
shuffle2    9,10,11,10
shuffle2    7,8,9,8
shuffle2    5,6,7,6
shuffle2    3,4,5,4
shuffle2    18,19,3,19


/* level 4 */
vmovdqa32	(_ZETAS_QINV+64)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+64)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20

butterfly	17,9,1,18,2,20
butterfly	16,8,1,18,2,20
butterfly	15,7,1,18,2,20
butterfly	14,6,1,18,2,20
butterfly	13,5,1,18,2,20
butterfly	12,4,1,18,2,20
butterfly	11,3,1,18,2,20
butterfly	10,19,1,18,2,20


/* level 5 */
vmovdqa32	(_ZETAS_QINV+80)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+80)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	17,13,1,18,2,20
butterfly	16,12,1,18,2,20
butterfly	15,11,1,18,2,20
butterfly	14,10,1,18,2,20

vmovdqa32	(_ZETAS_QINV+96)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+96)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	9,5,1,18,2,20
butterfly	8,4,1,18,2,20
butterfly	7,3,1,18,2,20
butterfly	6,19,1,18,2,20


/* level 6*/
vmovdqa32	(_ZETAS_QINV+112)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+112)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	17,15,1,18,2,20
butterfly	16,14,1,18,2,20

vmovdqa32	(_ZETAS_QINV+128)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+128)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	13,11,1,18,2,20
butterfly	12,10,1,18,2,20

vmovdqa32	(_ZETAS_QINV+144)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+144)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	9,7,1,18,2,20
butterfly	8,6,1,18,2,20

vmovdqa32	(_ZETAS_QINV+160)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+160)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	5,3,1,18,2,20
butterfly	4,19,1,18,2,20


/*level 7*/
vmovdqa32	(_ZETAS_QINV+176)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+176)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	17,16,1,18,2,20

vmovdqa32	(_ZETAS_QINV+192)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+192)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	15,14,1,18,2,20

vmovdqa32	(_ZETAS_QINV+208)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+208)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	13,12,1,18,2,20

vmovdqa32	(_ZETAS_QINV+224)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+224)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	11,10,1,18,2,20

vmovdqa32	(_ZETAS_QINV+240)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+240)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	9,8,1,18,2,20

vmovdqa32	(_ZETAS_QINV+256)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+256)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	7,6,1,18,2,20

vmovdqa32	(_ZETAS_QINV+272)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+272)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	5,4,1,18,2,20

vmovdqa32	(_ZETAS_QINV+288)*4(%rsi),%zmm1
vmovdqa32	(_ZETAS+288)*4(%rsi),%zmm2
vmovshdup   %zmm1,%zmm18
vmovshdup   %zmm2,%zmm20
butterfly	3,19,1,18,2,20

vmovdqa32		%zmm17,0+0(%rdi)
vmovdqa32		%zmm16,64+0(%rdi)
vmovdqa32		%zmm15,128+0(%rdi)
vmovdqa32	 	%zmm14,192+0(%rdi)
vmovdqa32		%zmm13,256+0(%rdi)
vmovdqa32		%zmm12,320+0(%rdi)
vmovdqa32		%zmm11,384+0(%rdi)
vmovdqa32	 	%zmm10,448+0(%rdi)
vmovdqa32		%zmm9,512+0(%rdi)
vmovdqa32		%zmm8,576+0(%rdi)
vmovdqa32		%zmm7,640+0(%rdi)
vmovdqa32	 	%zmm6,704+0(%rdi)
vmovdqa32		%zmm5,768+0(%rdi)
vmovdqa32		%zmm4,832+0(%rdi)
vmovdqa32		%zmm3,896+0(%rdi)
vmovdqa32	 	%zmm19,960+0(%rdi)

.endm

.text
.global cdecl(ntt_avx)
cdecl(ntt_avx):
vmovdqa32		_16XQ*4(%rsi),%zmm0

kmovw           (_ZETAS_QINV+5)*4(%rsi),%k1 
kmovw           (_ZETAS_QINV+6)*4(%rsi),%k2


levels0t7	

ret
