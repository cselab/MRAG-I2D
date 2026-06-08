#pragma once
#if defined(__aarch64__) || defined(__arm64__)
#define SSE2NEON_SUPPRESS_WARNINGS
#include <sse2neon.h>
#else
#include_next <pmmintrin.h>
#endif
