#pragma once


#if defined(__SSE2__) || _M_IX86_FP >= 2 || _M_X64 

#include "Simd4f_sse.hpp"

#else

#include "Simd4f_default.hpp"

#endif