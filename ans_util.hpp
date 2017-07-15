#pragma once

#include <cstdint>

namespace ans {

namespace constants {
    const uint32_t BLOCK_SIZE = 128;
    const uint8_t OUTPUT_BASE_LOG2 = 32;
    const uint64_t OUTPUT_BASE = 1ULL << OUTPUT_BASE_LOG2;
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

void flush_state(uint64_t state, uint8_t*& out, size_t num_bytes)
{
    for (size_t i = 0; i < num_bytes; i++) {
        uint8_t out_byte = state & 0xFF;
        out--;
        *out = out_byte;
        state >>= 8;
    }
}

uint64_t init_decoder(const uint8_t*& in, uint8_t num_bytes)
{
    uint64_t state = 0;
    for (size_t i = 0; i < num_bytes; i++) {
        uint8_t new_byte = *in++;
        state <<= 8;
        state = state + new_byte;
    }
    return state;
}

uint64_t next_power_of_two(uint64_t x)
{
    if (x == 0) {
        return 1;
    }
    uint32_t res = 63 - __builtin_clzll(x);
    return (1ULL << (res + 1));
}

bool is_power_of_two(uint64_t x) { return ((x != 0) && !(x & (x - 1))); }

template <uint8_t t_width>
void output_unit(uint8_t*& out, uint64_t& state)
{
    static_assert(t_width % 8 == 0, "can only write byte-multiple units");
    uint8_t w = t_width;
    while (w) {
        w -= 8;
        --out;
        *out = (uint8_t)(state & 0xFF);
        state = state >> 8;
    }
}

template <>
void output_unit<8>(uint8_t*& out, uint64_t& state)
{
    --out;
    *out = (uint8_t)(state & 0xFF);
    state = state >> 8;
}

template <>
void output_unit<16>(uint8_t*& out, uint64_t& state)
{
    out -= 2;
    uint16_t* out16 = reinterpret_cast<uint16_t*>(out);
    *out16 = (uint16_t)(state & 0xFFFF);
    state = state >> 16;
}

template <>
void output_unit<32>(uint8_t*& out, uint64_t& state)
{
    out -= 4;
    uint32_t* out32 = reinterpret_cast<uint32_t*>(out);
    *out32 = (uint32_t)(state & 0xFFFFFFFF);
    state = state >> 32;
}

template <uint8_t t_width>
void input_unit(const uint8_t*& in, uint64_t& state, std::size_t& enc_size)
{
    static_assert(t_width % 8 == 0, "can only read byte-multiple units");
    uint8_t w = t_width;
    while (w) {
        uint8_t new_byte = *in++;
        state = (state << t_width) | uint64_t(new_byte);
        w -= 8;
        enc_size--;
    }
}

template <>
void input_unit<8>(const uint8_t*& in, uint64_t& state, std::size_t& enc_size)
{
    uint8_t new_byte = *in++;
    state = (state << 8) | uint64_t(new_byte);
    enc_size--;
}

template <>
void input_unit<16>(const uint8_t*& in, uint64_t& state, std::size_t& enc_size)
{
    const uint16_t* in16 = reinterpret_cast<const uint16_t*>(in);
    uint64_t new_unit = *in16;
    state = (state << 16) | new_unit;
    in += 2;
    enc_size -= 2;
}

template <>
void input_unit<32>(const uint8_t*& in, uint64_t& state, std::size_t& enc_size)
{
    const uint32_t* in32 = reinterpret_cast<const uint32_t*>(in);
    uint64_t new_unit = *in32;
    state = (state << 32) | new_unit;
    in += 4;
    enc_size -= 4;
}

uint8_t state_bytes(uint64_t state)
{
    return 8 - (__builtin_clzll(state) >> 3);
}

uint8_t pack_two_4bit_nums(uint8_t a, uint8_t b)
{
    return (a << 4) + b;
}

std::pair<uint8_t, uint8_t> unpack_two_4bit_nums(uint8_t x)
{
    return { (x >> 4), (x & 15) };
}
}