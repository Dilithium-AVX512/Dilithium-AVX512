#include "consts.h"
.include "shuffle.inc"

.macro butterfly l,h,zl0=1,zl1=1,zh0=2,zh1=2
vpsubd		%zmm\l,%zmm\h,%zmm28
vpaddd		%zmm\h,%zmm\l,%zmm\l

vpmuldq		%zmm\zl0,%zmm28,%zmm29
vmovshdup	%zmm28,%zmm\h
vpmuldq		%zmm\zl1,%zmm\h,%zmm30

vpmuldq		%zmm\zh0,%zmm28,%zmm28
vpmuldq		%zmm\zh1,%zmm\h,%zmm\h

vpmuldq		%zmm0,%zmm29,%zmm29
vpmuldq		%zmm0,%zmm30,%zmm30

vpsubd		%zmm29,%zmm28,%zmm28
vpsubd		%zmm30,%zmm\h,%zmm\h

vmovshdup	%zmm28,%zmm28
vpblendmd	%zmm\h,%zmm28,%zmm\h {%k1}
.endm

.macro levels0t7 
vmovdqa32		0+0(%rdi),%zmm4
vmovdqa32		64+0(%rdi),%zmm5
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


/* level 0 */
vpermd		(_ZETAS_QINV+304-16)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+304-16)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	4,5,1,3,2,20

vpermd		(_ZETAS_QINV+304-16*2)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+304-16*2)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	6,7,1,3,2,20

vpermd		(_ZETAS_QINV+304-16*3)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+304-16*3)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	8,9,1,3,2,20

vpermd		(_ZETAS_QINV+304-16*4)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+304-16*4)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	10,11,1,3,2,20

vpermd		(_ZETAS_QINV+304-16*5)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+304-16*5)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	12,13,1,3,2,20

vpermd		(_ZETAS_QINV+304-16*6)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+304-16*6)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	14,15,1,3,2,20

vpermd		(_ZETAS_QINV+304-16*7)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+304-16*7)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	16,17,1,3,2,20

vpermd		(_ZETAS_QINV+304-16*8)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+304-16*8)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	18,19,1,3,2,20


/*level 1*/
vpermd		(_ZETAS_QINV+176-16)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+176-16)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	4,6,1,3,2,20
butterfly	5,7,1,3,2,20

vpermd		(_ZETAS_QINV+176-16*2)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+176-16*2)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	8,10,1,3,2,20
butterfly	9,11,1,3,2,20

vpermd		(_ZETAS_QINV+176-16*3)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+176-16*3)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	12,14,1,3,2,20
butterfly	13,15,1,3,2,20

vpermd		(_ZETAS_QINV+176-16*4)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+176-16*4)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	16,18,1,3,2,20
butterfly	17,19,1,3,2,20


/* level 2 */
vpermd		(_ZETAS_QINV+112-16)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+112-16)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	4,8,1,3,2,20
butterfly	5,9,1,3,2,20
butterfly	6,10,1,3,2,20
butterfly	7,11,1,3,2,20

vpermd		(_ZETAS_QINV+112-16*2)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+112-16*2)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	12,16,1,3,2,20
butterfly	13,17,1,3,2,20
butterfly	14,18,1,3,2,20
butterfly	15,19,1,3,2,20


/* level 3 */
vpermd		(_ZETAS_QINV+64)*4(%rsi),%zmm21,%zmm3
vpermd		(_ZETAS+64)*4(%rsi),%zmm21,%zmm20
vmovshdup	%zmm3,%zmm1
vmovshdup	%zmm20,%zmm2

butterfly	4,12,1,3,2,20
butterfly	5,13,1,3,2,20
butterfly	6,14,1,3,2,20
butterfly	7,15,1,3,2,20
butterfly	8,16,1,3,2,20
butterfly	9,17,1,3,2,20
butterfly	10,18,1,3,2,20
butterfly	11,19,1,3,2,20


/* level 4 */
shuffle2   4,5,3,5
shuffle2   6,7,4,7
shuffle2   8,9,6,9
shuffle2   10,11,8,11
shuffle2   12,13,10,13
shuffle2   14,15,12,15
shuffle2   16,17,14,17
shuffle2   18,19,16,19

vpermd		(_ZETAS_QINV+48)*4(%rsi),%zmm21,%zmm1
vpermd		(_ZETAS+48)*4(%rsi),%zmm21,%zmm2

butterfly	3,5
butterfly	4,7
butterfly	6,9
butterfly	8,11
butterfly	10,13
butterfly	12,15
butterfly	14,17
butterfly	16,19


/* level 5 */
shuffle4   3,4,18,4
shuffle4   6,8,3,8
shuffle4   10,12,6,12
shuffle4   14,16,10,16
shuffle4   5,7,14,7
shuffle4   9,11,5,11
shuffle4   13,15,9,15
shuffle4   17,19,13,19

vpermd		(_ZETAS_QINV+32)*4(%rsi),%zmm21,%zmm1
vpermd		(_ZETAS+32)*4(%rsi),%zmm21,%zmm2

butterfly	18,4
butterfly	3,8
butterfly	6,12
butterfly	10,16
butterfly	14,7
butterfly	5,11
butterfly	9,15
butterfly	13,19


/*level 6*/
shuffle8   18,3,17,3
shuffle8   6,10,18,10
shuffle8   14,5,6,5
shuffle8   9,13,14,13
shuffle8   4,8,9,8
shuffle8   12,16,4,16
shuffle8   7,11,12,11
shuffle8   15,19,7,19

vpermd		(_ZETAS_QINV+16)*4(%rsi),%zmm21,%zmm1
vpermd		(_ZETAS+16)*4(%rsi),%zmm21,%zmm2

butterfly   17,3
butterfly	18,10
butterfly	6,5
butterfly	14,13
butterfly	9,8
butterfly	4,16
butterfly	12,11
butterfly	7,19


/*level 7*/
shuffle16   17,18,15,18
shuffle16   6,14,17,14
shuffle16   9,4,6,4
shuffle16   12,7,9,7
shuffle16   3,10,12,10
shuffle16   5,13,3,13
shuffle16   8,16,5,16
shuffle16   11,19,8,19

vpbroadcastd	(_ZETAS_QINV+0)*4(%rsi),%zmm1
vpbroadcastd	(_ZETAS+0)*4(%rsi),%zmm2

butterfly   15,18
butterfly	17,14
butterfly	6,4
butterfly	9,7
butterfly	12,10
butterfly	3,13
butterfly	5,16
butterfly	8,19

vmovdqa32		%zmm18,512+0(%rdi)
vmovdqa32		%zmm14,576+0(%rdi)
vmovdqa32		%zmm4,640+0(%rdi)
vmovdqa32	 	%zmm7,704+0(%rdi)
vmovdqa32		%zmm10,768+0(%rdi)
vmovdqa32		%zmm13,832+0(%rdi)
vmovdqa32		%zmm16,896+0(%rdi)
vmovdqa32	 	%zmm19,960+0(%rdi)

vmovdqa32		(_16XDIV_QINV)*4(%rsi),%zmm1
vmovdqa32		(_16XDIV)*4(%rsi),%zmm2
vpmuldq		    %zmm1,%zmm15,%zmm22
vpmuldq		    %zmm1,%zmm17,%zmm23
vpmuldq		    %zmm1,%zmm6,%zmm24
vpmuldq		    %zmm1,%zmm9,%zmm25
vmovshdup	    %zmm15,%zmm18
vmovshdup	    %zmm17,%zmm14
vmovshdup	    %zmm6,%zmm4
vmovshdup	    %zmm9,%zmm7
vpmuldq		    %zmm1,%zmm18,%zmm26
vpmuldq		    %zmm1,%zmm14,%zmm27
vpmuldq		    %zmm1,%zmm4,%zmm28
vpmuldq		    %zmm1,%zmm7,%zmm29
vpmuldq		    %zmm2,%zmm15,%zmm15
vpmuldq		    %zmm2,%zmm17,%zmm17
vpmuldq		    %zmm2,%zmm6,%zmm6
vpmuldq		    %zmm2,%zmm9,%zmm9
vpmuldq		    %zmm2,%zmm18,%zmm18
vpmuldq		    %zmm2,%zmm14,%zmm14
vpmuldq		    %zmm2,%zmm4,%zmm4
vpmuldq		    %zmm2,%zmm7,%zmm7
vpmuldq		    %zmm0,%zmm22,%zmm22
vpmuldq		    %zmm0,%zmm23,%zmm23
vpmuldq		    %zmm0,%zmm24,%zmm24
vpmuldq		    %zmm0,%zmm25,%zmm25
vpmuldq		    %zmm0,%zmm26,%zmm26
vpmuldq		    %zmm0,%zmm27,%zmm27
vpmuldq		    %zmm0,%zmm28,%zmm28
vpmuldq		    %zmm0,%zmm29,%zmm29
vpsubd		    %zmm22,%zmm15,%zmm15
vpsubd		    %zmm23,%zmm17,%zmm17
vpsubd		    %zmm24,%zmm6,%zmm6
vpsubd		    %zmm25,%zmm9,%zmm9
vpsubd		    %zmm26,%zmm18,%zmm18
vpsubd		    %zmm27,%zmm14,%zmm14
vpsubd		    %zmm28,%zmm4,%zmm4
vpsubd		    %zmm29,%zmm7,%zmm7
vmovshdup	    %zmm15,%zmm15
vmovshdup	    %zmm17,%zmm17
vmovshdup	    %zmm6,%zmm6
vmovshdup	    %zmm9,%zmm9
vpblendmd	    %zmm18,%zmm15,%zmm15 {%k1}
vpblendmd	    %zmm14,%zmm17,%zmm17 {%k1}
vpblendmd	    %zmm4,%zmm6,%zmm6 {%k1}
vpblendmd	    %zmm7,%zmm9,%zmm9 {%k1}



vpmuldq		    %zmm1,%zmm12,%zmm22
vpmuldq		    %zmm1,%zmm3,%zmm23
vpmuldq		    %zmm1,%zmm5,%zmm24
vpmuldq		    %zmm1,%zmm8,%zmm25
vmovshdup	    %zmm12,%zmm18
vmovshdup	    %zmm3,%zmm14
vmovshdup	    %zmm5,%zmm4
vmovshdup	    %zmm8,%zmm7
vpmuldq		%zmm1,%zmm18,%zmm26
vpmuldq		%zmm1,%zmm14,%zmm27
vpmuldq		%zmm1,%zmm4,%zmm28
vpmuldq		%zmm1,%zmm7,%zmm29
vpmuldq		%zmm2,%zmm12,%zmm12
vpmuldq		%zmm2,%zmm3,%zmm3
vpmuldq		%zmm2,%zmm5,%zmm5
vpmuldq		%zmm2,%zmm8,%zmm8
vpmuldq		%zmm2,%zmm18,%zmm18
vpmuldq		%zmm2,%zmm14,%zmm14
vpmuldq		%zmm2,%zmm4,%zmm4
vpmuldq		%zmm2,%zmm7,%zmm7
vpmuldq		%zmm0,%zmm22,%zmm22
vpmuldq		%zmm0,%zmm23,%zmm23
vpmuldq		%zmm0,%zmm24,%zmm24
vpmuldq		%zmm0,%zmm25,%zmm25
vpmuldq		%zmm0,%zmm26,%zmm26
vpmuldq		%zmm0,%zmm27,%zmm27
vpmuldq		%zmm0,%zmm28,%zmm28
vpmuldq		%zmm0,%zmm29,%zmm29
vpsubd		%zmm22,%zmm12,%zmm12
vpsubd		%zmm23,%zmm3,%zmm3
vpsubd		%zmm24,%zmm5,%zmm5
vpsubd		%zmm25,%zmm8,%zmm8
vpsubd		%zmm26,%zmm18,%zmm18
vpsubd		%zmm27,%zmm14,%zmm14
vpsubd		%zmm28,%zmm4,%zmm4
vpsubd		%zmm29,%zmm7,%zmm7
vmovshdup	%zmm12,%zmm12
vmovshdup	%zmm3,%zmm3
vmovshdup	%zmm5,%zmm5
vmovshdup	%zmm8,%zmm8
vpblendmd	%zmm18,%zmm12,%zmm12 {%k1}
vpblendmd	%zmm14,%zmm3,%zmm3 {%k1}
vpblendmd	%zmm4,%zmm5,%zmm5 {%k1}
vpblendmd	%zmm7,%zmm8,%zmm8 {%k1}

vmovdqa32		%zmm15,0+0(%rdi)
vmovdqa32		%zmm17,64+0(%rdi)
vmovdqa32		%zmm6,128+0(%rdi)
vmovdqa32	 	%zmm9,192+0(%rdi)
vmovdqa32		%zmm12,256+0(%rdi)
vmovdqa32		%zmm3,320+0(%rdi)
vmovdqa32		%zmm5,384+0(%rdi)
vmovdqa32	 	%zmm8,448+0(%rdi)

.endm


.text
.global cdecl(invntt_avx)
cdecl(invntt_avx):
vmovdqa32		_16XQ*4(%rsi),%zmm0

vmovdqa32	    _zmm21*4(%rsi),%zmm21


kmovw           (_ZETAS_QINV+5)*4(%rsi),%k1 


kmovw           (_ZETAS_QINV+6)*4(%rsi),%k2

levels0t7	

ret
