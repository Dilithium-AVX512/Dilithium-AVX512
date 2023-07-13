#include "params.h"
#include "consts.h"

.text
.global cdecl(pointwise_avx)
cdecl(pointwise_avx):
#consts
vmovdqa32	_16XQINV*4(%rcx),%zmm0
vmovdqa32	_16XQ*4(%rcx),%zmm1

kmovw           (_ZETAS_QINV+5)*4(%rcx),%k1 

xor		%eax,%eax
_looptop1:
#load
vmovdqa32		(%rsi),%zmm2
vmovdqa32		64(%rsi),%zmm4
vmovdqa32		128(%rsi),%zmm6
vmovdqa32		192(%rsi),%zmm8
vmovdqa32		256(%rsi),%zmm10
vmovdqa32		320(%rsi),%zmm12
vmovdqa32		384(%rsi),%zmm14

vmovdqa32		(%rdx),%zmm18
vmovdqa32		64(%rdx),%zmm20
vmovdqa32		128(%rdx),%zmm22
vmovdqa32		192(%rdx),%zmm24
vmovdqa32		256(%rdx),%zmm26
vmovdqa32		320(%rdx),%zmm28
vmovdqa32		384(%rdx),%zmm30

vmovshdup	%zmm2,%zmm3
vmovshdup	%zmm4,%zmm5
vmovshdup	%zmm6,%zmm7
vmovshdup	%zmm8,%zmm9
vmovshdup	%zmm10,%zmm11
vmovshdup	%zmm12,%zmm13
vmovshdup	%zmm14,%zmm15
vmovshdup	%zmm18,%zmm19
vmovshdup	%zmm20,%zmm21
vmovshdup	%zmm22,%zmm23
vmovshdup	%zmm24,%zmm25
vmovshdup	%zmm26,%zmm27
vmovshdup	%zmm28,%zmm29
vmovshdup	%zmm30,%zmm31

#mul
vpmuldq		%zmm2,%zmm18,%zmm2
vpmuldq		%zmm3,%zmm19,%zmm3
vpmuldq		%zmm4,%zmm20,%zmm4
vpmuldq		%zmm5,%zmm21,%zmm5
vpmuldq		%zmm6,%zmm22,%zmm6
vpmuldq		%zmm7,%zmm23,%zmm7
vpmuldq		%zmm8,%zmm24,%zmm8
vpmuldq		%zmm9,%zmm25,%zmm9
vpmuldq		%zmm10,%zmm26,%zmm10
vpmuldq		%zmm11,%zmm27,%zmm11
vpmuldq		%zmm12,%zmm28,%zmm12
vpmuldq		%zmm13,%zmm29,%zmm13
vpmuldq		%zmm14,%zmm30,%zmm14
vpmuldq		%zmm15,%zmm31,%zmm15

#reduce
vpmuldq		%zmm0,%zmm2,%zmm18
vpmuldq		%zmm0,%zmm3,%zmm19
vpmuldq		%zmm0,%zmm4,%zmm20
vpmuldq		%zmm0,%zmm5,%zmm21
vpmuldq		%zmm0,%zmm6,%zmm22
vpmuldq		%zmm0,%zmm7,%zmm23
vpmuldq		%zmm0,%zmm8,%zmm24
vpmuldq		%zmm0,%zmm9,%zmm25
vpmuldq		%zmm0,%zmm10,%zmm26
vpmuldq		%zmm0,%zmm11,%zmm27
vpmuldq		%zmm0,%zmm12,%zmm28
vpmuldq		%zmm0,%zmm13,%zmm29
vpmuldq		%zmm0,%zmm14,%zmm30
vpmuldq		%zmm0,%zmm15,%zmm31

vpmuldq		%zmm1,%zmm18,%zmm18
vpmuldq		%zmm1,%zmm19,%zmm19
vpmuldq		%zmm1,%zmm20,%zmm20
vpmuldq		%zmm1,%zmm21,%zmm21
vpmuldq		%zmm1,%zmm22,%zmm22
vpmuldq		%zmm1,%zmm23,%zmm23
vpmuldq		%zmm1,%zmm24,%zmm24
vpmuldq		%zmm1,%zmm25,%zmm25
vpmuldq		%zmm1,%zmm26,%zmm26
vpmuldq		%zmm1,%zmm27,%zmm27
vpmuldq		%zmm1,%zmm28,%zmm28
vpmuldq		%zmm1,%zmm29,%zmm29
vpmuldq		%zmm1,%zmm30,%zmm30
vpmuldq		%zmm1,%zmm31,%zmm31

vpsubq		%zmm18,%zmm2,%zmm2
vpsubq		%zmm19,%zmm3,%zmm3
vpsubq		%zmm20,%zmm4,%zmm4
vpsubq		%zmm21,%zmm5,%zmm5
vpsubq		%zmm22,%zmm6,%zmm6
vpsubq		%zmm23,%zmm7,%zmm7
vpsubq		%zmm24,%zmm8,%zmm8
vpsubq		%zmm25,%zmm9,%zmm9
vpsubq		%zmm26,%zmm10,%zmm10
vpsubq		%zmm27,%zmm11,%zmm11
vpsubq		%zmm28,%zmm12,%zmm12
vpsubq		%zmm29,%zmm13,%zmm13
vpsubq		%zmm30,%zmm14,%zmm14
vpsubq		%zmm31,%zmm15,%zmm15


vmovshdup	%zmm2,%zmm2
vmovshdup	%zmm4,%zmm4
vmovshdup	%zmm6,%zmm6
vmovshdup	%zmm8,%zmm8
vmovshdup	%zmm10,%zmm10
vmovshdup	%zmm12,%zmm12
vmovshdup	%zmm14,%zmm14

#store
vpblendmd	%zmm3,%zmm2,%zmm2 {%k1}
vpblendmd	%zmm5,%zmm4,%zmm4 {%k1}
vpblendmd	%zmm7,%zmm6,%zmm6 {%k1}
vpblendmd	%zmm9,%zmm8,%zmm8 {%k1}
vpblendmd	%zmm11,%zmm10,%zmm10 {%k1}
vpblendmd	%zmm13,%zmm12,%zmm12 {%k1}
vpblendmd	%zmm15,%zmm14,%zmm14 {%k1}

vmovdqa32		%zmm2,(%rdi)
vmovdqa32		%zmm4,64(%rdi)
vmovdqa32		%zmm6,128(%rdi)
vmovdqa32		%zmm8,192(%rdi)
vmovdqa32		%zmm10,256(%rdi)
vmovdqa32		%zmm12,320(%rdi)
vmovdqa32		%zmm14,384(%rdi)

add		$448,%rdi
add		$448,%rsi
add		$448,%rdx
add		$1,%eax
cmp		$2,%eax
jb 		_looptop1

vmovdqa32	(%rsi),%zmm2
vmovdqa32	64(%rsi),%zmm4
vmovdqa32	(%rdx),%zmm10
vmovdqa32	64(%rdx),%zmm12
vmovshdup	%zmm2,%zmm3
vmovshdup	%zmm4,%zmm5
vmovshdup	%zmm10,%zmm11
vmovshdup	%zmm12,%zmm13

#mul
vpmuldq		%zmm2,%zmm10,%zmm2
vpmuldq		%zmm3,%zmm11,%zmm3
vpmuldq		%zmm4,%zmm12,%zmm4
vpmuldq		%zmm5,%zmm13,%zmm5

#reduce
vpmuldq		%zmm0,%zmm2,%zmm10
vpmuldq		%zmm0,%zmm3,%zmm11
vpmuldq		%zmm0,%zmm4,%zmm12
vpmuldq		%zmm0,%zmm5,%zmm13

vpmuldq		%zmm1,%zmm10,%zmm10
vpmuldq		%zmm1,%zmm11,%zmm11
vpmuldq		%zmm1,%zmm12,%zmm12
vpmuldq		%zmm1,%zmm13,%zmm13

vpsubq		%zmm10,%zmm2,%zmm2
vpsubq		%zmm11,%zmm3,%zmm3
vpsubq		%zmm12,%zmm4,%zmm4
vpsubq		%zmm13,%zmm5,%zmm5

vmovshdup	%zmm2,%zmm2
vmovshdup	%zmm4,%zmm4

#store
vpblendmd	%zmm3,%zmm2,%zmm2 {%k1}
vpblendmd	%zmm5,%zmm4,%zmm4 {%k1}
vmovdqa32   %zmm2,(%rdi)
vmovdqa32	%zmm4,64(%rdi)

ret




.macro pointwise off
#load
vmovdqa32	\off(%rsi),%zmm10
vmovdqa32	\off+64(%rsi),%zmm12
vmovdqa32	\off+128(%rsi),%zmm14
vmovdqa32	\off+192(%rsi),%zmm16
vmovdqa32	\off(%rdx),%zmm18
vmovdqa32	\off+64(%rdx),%zmm20
vmovdqa32	\off+128(%rdx),%zmm22
vmovdqa32	\off+192(%rdx),%zmm24
vmovshdup	%zmm10,%zmm11
vmovshdup	%zmm12,%zmm13
vmovshdup	%zmm14,%zmm15
vmovshdup	%zmm16,%zmm17
vmovshdup	%zmm18,%zmm19
vmovshdup	%zmm20,%zmm21
vmovshdup	%zmm22,%zmm23
vmovshdup	%zmm24,%zmm25

#mul
vpmuldq		%zmm10,%zmm18,%zmm10
vpmuldq		%zmm11,%zmm19,%zmm11
vpmuldq		%zmm12,%zmm20,%zmm12
vpmuldq		%zmm13,%zmm21,%zmm13
vpmuldq		%zmm14,%zmm22,%zmm14
vpmuldq		%zmm15,%zmm23,%zmm15
vpmuldq		%zmm16,%zmm24,%zmm16
vpmuldq		%zmm17,%zmm25,%zmm17
.endm

.macro acc
vpaddq		%zmm10,%zmm2,%zmm2
vpaddq		%zmm11,%zmm3,%zmm3
vpaddq		%zmm12,%zmm4,%zmm4
vpaddq		%zmm13,%zmm5,%zmm5
vpaddq		%zmm14,%zmm6,%zmm6
vpaddq		%zmm15,%zmm7,%zmm7
vpaddq		%zmm16,%zmm8,%zmm8
vpaddq		%zmm17,%zmm9,%zmm9
.endm

.global cdecl(pointwise_acc_avx)
cdecl(pointwise_acc_avx):
#consts
vmovdqa32		_16XQINV*4(%rcx),%zmm0
vmovdqa32		_16XQ*4(%rcx),%zmm1
kmovw           (_ZETAS_QINV+5)*4(%rcx),%k1 

xor		%eax,%eax
_looptop2:
pointwise	0

#mov
vmovdqa32		%zmm10,%zmm2
vmovdqa32		%zmm11,%zmm3
vmovdqa32		%zmm12,%zmm4
vmovdqa32		%zmm13,%zmm5
vmovdqa32		%zmm14,%zmm6
vmovdqa32		%zmm15,%zmm7
vmovdqa32		%zmm16,%zmm8
vmovdqa32		%zmm17,%zmm9

pointwise	1024
acc

#if L >= 3
pointwise	2048
acc
#endif

#if L >= 4
pointwise	3072
acc
#endif

#if L >= 5
pointwise	4096
acc
#endif

#if L >= 6
pointwise	5120
acc
#endif

#if L >= 7
pointwise	6144
acc
#endif

#reduce
vpmuldq		%zmm0,%zmm2,%zmm10
vpmuldq		%zmm0,%zmm3,%zmm11
vpmuldq		%zmm0,%zmm4,%zmm12
vpmuldq		%zmm0,%zmm5,%zmm13
vpmuldq		%zmm0,%zmm6,%zmm14
vpmuldq		%zmm0,%zmm7,%zmm15
vpmuldq		%zmm0,%zmm8,%zmm16
vpmuldq		%zmm0,%zmm9,%zmm17

vpmuldq		%zmm1,%zmm10,%zmm10
vpmuldq		%zmm1,%zmm11,%zmm11
vpmuldq		%zmm1,%zmm12,%zmm12
vpmuldq		%zmm1,%zmm13,%zmm13
vpmuldq		%zmm1,%zmm14,%zmm14
vpmuldq		%zmm1,%zmm15,%zmm15
vpmuldq		%zmm1,%zmm16,%zmm16
vpmuldq		%zmm1,%zmm17,%zmm17

vpsubq		%zmm10,%zmm2,%zmm2
vpsubq		%zmm11,%zmm3,%zmm3
vpsubq		%zmm12,%zmm4,%zmm4
vpsubq		%zmm13,%zmm5,%zmm5
vpsubq		%zmm14,%zmm6,%zmm6
vpsubq		%zmm15,%zmm7,%zmm7
vpsubq		%zmm16,%zmm8,%zmm8
vpsubq		%zmm17,%zmm9,%zmm9

vmovshdup	%zmm2,%zmm2
vmovshdup	%zmm4,%zmm4
vmovshdup	%zmm6,%zmm6
vmovshdup	%zmm8,%zmm8

#store
vpblendmd	%zmm3,%zmm2,%zmm2 {%k1}
vpblendmd	%zmm5,%zmm4,%zmm4 {%k1}
vpblendmd	%zmm7,%zmm6,%zmm6 {%k1}
vpblendmd	%zmm9,%zmm8,%zmm8 {%k1}

vmovdqa32	%zmm2,(%rdi)
vmovdqa32	%zmm4,64(%rdi)
vmovdqa32	%zmm6,128(%rdi)
vmovdqa32	%zmm8,192(%rdi)

add		$256,%rsi
add		$256,%rdx
add		$256,%rdi
add		$1,%eax
cmp		$4,%eax
jb _looptop2

ret

