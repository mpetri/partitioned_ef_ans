#pragma once

#include <cstdint>

namespace ans {

namespace constants {
    const uint32_t BLOCK_SIZE = 128;
}

uint8_t magnitude(uint32_t x)
{
    uint64_t y = x;
    if (x == 1)
        return 0;
    uint32_t res = 63 - __builtin_clzll(y);
    if ((1ULL << res) == y)
        return res;
    return res + 1;
}
}