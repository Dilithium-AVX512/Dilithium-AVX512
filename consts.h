#ifndef CONSTS_H
#define CONSTS_H

#include "params.h"

#define _16XQ          0
#define _16XQINV       16
#define _16XDIV_QINV   32
#define _16XDIV        48
#define _zmm21        64
#define _ZETAS_QINV   80
#define _ZETAS        384


/* The C ABI on MacOS exports all symbols with a leading
 * underscore. This means that any symbols we refer to from
 * C files (functions) can't be found, and all symbols we
 * refer to from ASM also can't be found.
 *
 * This define helps us get around this
 */
#if defined(__WIN32__) || defined(__APPLE__)
#define decorate(s) _##s
#define _cdecl(s) decorate(s)
#define cdecl(s) _cdecl(DILITHIUM_NAMESPACE(##s))
#else
#define cdecl(s) DILITHIUM_NAMESPACE(##s)
#endif

#ifndef __ASSEMBLER__

#include "align.h"

typedef ALIGNED_INT32(688) qdata_t;

#define qdata DILITHIUM_NAMESPACE(qdata)
extern const qdata_t qdata;

#endif
#endif
