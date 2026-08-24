#pragma once

#include "zk_config.h"

// backend implementation
#if defined(AUTOHOG_ZK_BACKEND_mcl)
#include "mcl_wrapper.h"

#elif defined(AUTOHOG_ZK_BACKEND_blst)

// In C++17, std::byte conflicts with blst's byte typedef.
#define byte blst_byte_t
#include <blst.h>
#undef byte

#else
#error "Unknown ZK backend selection"
#endif
