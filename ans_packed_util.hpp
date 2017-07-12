#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <utility>

namespace ans_packed {

namespace constants {
    const uint32_t BLOCK_SIZE = 128;
    const uint8_t MAX_MAG = 32;
    const uint8_t NUM_MODELS = 16;
    const std::array<uint8_t, NUM_MODELS> SEL2MAG{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12,
        14, 16, 19, 22, 32 };
    const std::array<uint8_t, MAX_MAG + 1> MAG2SEL{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
        10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
    const uint64_t TOPFREQ = 1048576;
    const uint8_t OUTPUT_BASE_LOG2 = 32;
    const uint64_t OUTPUT_BASE = 1ULL << OUTPUT_BASE_LOG2;
    const uint64_t NORM_LOWER_BOUND = 1ULL << 24;
}

struct mag_enc_table_entry {
    uint32_t freq;
    uint64_t base;
    uint64_t SUB;
};

std::ostream& operator<<(std::ostream& os, const mag_enc_table_entry& t)
{
    os << "<f=" << t.freq << ",b=" << t.base << ",SUB=" << t.SUB << ">";
    return os;
}

struct mag_dec_table_entry {
    uint32_t freq;
    uint64_t offset;
    uint32_t sym;
};

struct enc_model {
    uint64_t M = 0; // frame size
    uint8_t log2_M = 0;
    uint64_t mask_M = 0;
    uint64_t norm_lower_bound = 0;
    uint32_t max_value = 0;
    mag_enc_table_entry table[0];
};

void print_enc_model(const enc_model* m)
{
    std::cout << "ENC_MODEL M=" << m->M
              << " log2_M = " << (int)m->log2_M
              << " mask_M = " << m->mask_M
              << " norm_lower_bound = " << m->norm_lower_bound
              << " max_value = " << (int)m->max_value << "\n";
    std::cout << "TABLE = [";
    for (size_t i = 0; i < m->max_value; i++) {
        std::cout << m->table[i];
    }
    std::cout << "]\n";
}

struct mag_table {
    uint32_t max_value;
    uint64_t counts[constants::MAX_MAG + 1];
};

void print_mag_table(const mag_table* tb, std::string name)
{
    std::cout << name << " max_value = " << tb->max_value << " COUNTS = (";
    for (size_t i = 0; i <= constants::MAX_MAG; i++) {
        std::cout << "<" << i << "," << tb->counts[i] << ">";
    }
    std::cout << ")" << std::endl;
}

struct dec_model {
    uint64_t M = 0; // frame size
    uint8_t log2_M = 0;
    uint64_t mask_M = 0;
    uint64_t norm_lower_bound = 0;
    mag_dec_table_entry table[0];
};

uint8_t state_bytes(uint64_t state)
{
    if (state < (1ULL << 8))
        return 1;
    if (state < (1ULL << 16))
        return 2;
    if (state < (1ULL << 24))
        return 3;
    if (state < (1ULL << 32))
        return 4;
    if (state < (1ULL << 40))
        return 5;
    if (state < (1ULL << 48))
        return 6;
    if (state < (1ULL << 56))
        return 7;
    return 8;
}

uint8_t pack_two_4bit_nums(uint8_t a, uint8_t b)
{
    return (a << 4) + b;
}

std::pair<uint8_t, uint8_t> unpack_two_4bit_nums(uint8_t x)
{
    return { (x >> 4), (x & 15) };
}

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

inline uint8_t vb_size(uint64_t x)
{
    if (x < (1ULL << 7)) {
        return 1;
    } else if (x < (1ULL << 14)) {
        return 2;
    } else if (x < (1ULL << 21)) {
        return 3;
    } else if (x < (1ULL << 28)) {
        return 4;
    } else if (x < (1ULL << 35)) {
        return 5;
    } else if (x < (1ULL << 42)) {
        return 6;
    } else if (x < (1ULL << 49)) {
        return 7;
    } else if (x < (1ULL << 56)) {
        return 8;
    }
    return 9;
}

inline uint64_t vbyte_decode_u64(const uint8_t*& input)
{
    uint64_t x = 0;
    uint64_t shift = 0;
    while (true) {
        uint8_t c = *input++;
        x += (uint64_t(c & 127) << shift);
        if (!(c & 128)) {
            return x;
        }
        shift += 7;
    }
    return x;
}

template <uint32_t i>
inline uint8_t extract7bits(const uint64_t val)
{
    uint8_t v = static_cast<uint8_t>((val >> (7 * i)) & ((1ULL << 7) - 1));
    return v;
}

template <uint32_t i>
inline uint8_t extract7bitsmaskless(const uint64_t val)
{
    uint8_t v = static_cast<uint8_t>((val >> (7 * i)));
    return v;
}

inline void vbyte_encode_u64(uint8_t*& out, uint64_t x)
{
    if (x < (1ULL << 7)) {
        *out++ = static_cast<uint8_t>(x & 127);
    } else if (x < (1ULL << 14)) {
        *out++ = extract7bits<0>(x) | 128;
        *out++ = extract7bitsmaskless<1>(x) & 127;
    } else if (x < (1ULL << 21)) {
        *out++ = extract7bits<0>(x) | 128;
        *out++ = extract7bits<1>(x) | 128;
        *out++ = extract7bitsmaskless<2>(x) & 127;
    } else if (x < (1ULL << 28)) {
        *out++ = extract7bits<0>(x) | 128;
        *out++ = extract7bits<1>(x) | 128;
        *out++ = extract7bits<2>(x) | 128;
        *out++ = extract7bitsmaskless<3>(x) & 127;
    } else if (x < (1ULL << 35)) {
        *out++ = extract7bits<0>(x) | 128;
        *out++ = extract7bits<1>(x) | 128;
        *out++ = extract7bits<2>(x) | 128;
        *out++ = extract7bits<3>(x) | 128;
        *out++ = extract7bitsmaskless<4>(x) & 127;
    } else if (x < (1ULL << 42)) {
        *out++ = extract7bits<0>(x) | 128;
        *out++ = extract7bits<1>(x) | 128;
        *out++ = extract7bits<2>(x) | 128;
        *out++ = extract7bits<3>(x) | 128;
        *out++ = extract7bits<4>(x) | 128;
        *out++ = extract7bitsmaskless<5>(x) & 127;
    } else if (x < (1ULL << 49)) {
        *out++ = extract7bits<0>(x) | 128;
        *out++ = extract7bits<1>(x) | 128;
        *out++ = extract7bits<2>(x) | 128;
        *out++ = extract7bits<3>(x) | 128;
        *out++ = extract7bits<4>(x) | 128;
        *out++ = extract7bits<5>(x) | 128;
        *out++ = extract7bitsmaskless<6>(x) & 127;
    } else if (x < (1ULL << 56)) {
        *out++ = extract7bits<0>(x) | 128;
        *out++ = extract7bits<1>(x) | 128;
        *out++ = extract7bits<2>(x) | 128;
        *out++ = extract7bits<3>(x) | 128;
        *out++ = extract7bits<4>(x) | 128;
        *out++ = extract7bits<5>(x) | 128;
        *out++ = extract7bits<6>(x) | 128;
        *out++ = extract7bitsmaskless<7>(x) & 127;
    } else {
        *out++ = extract7bits<0>(x) | 128;
        *out++ = extract7bits<1>(x) | 128;
        *out++ = extract7bits<2>(x) | 128;
        *out++ = extract7bits<3>(x) | 128;
        *out++ = extract7bits<4>(x) | 128;
        *out++ = extract7bits<5>(x) | 128;
        *out++ = extract7bits<6>(x) | 128;
        *out++ = extract7bits<7>(x) | 128;
        *out++ = extract7bitsmaskless<8>(x) & 127;
    }
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

uint32_t max_val_in_mag(uint8_t mag, uint32_t max_val = 0)
{
    uint32_t maxv = 1;
    if (mag != 0)
        maxv = (1ULL << (mag));
    if (maxv > max_val)
        maxv = max_val;
    return maxv;
}

uint32_t min_val_in_mag(uint8_t mag)
{
    if (mag == 0)
        return 1;
    return (1ULL << (mag - 1)) + 1;
}

uint32_t uniq_vals_in_mag(uint8_t mag, uint32_t max_val = 0)
{
    return max_val_in_mag(mag, max_val) - min_val_in_mag(mag) + 1;
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

static mag_table* normalize_counts(const mag_table* table)
{
    print_mag_table(table, "initial_freqs");
    mag_table* nfreqs = new mag_table;
    nfreqs->max_value = table->max_value;
    uint64_t initial_sum = 0;
    for (size_t i = 0; i <= constants::MAX_MAG; i++) {
        initial_sum += table->counts[i] * uniq_vals_in_mag(i, table->max_value);
        nfreqs->counts[i] = table->counts[i];
    }

    uint8_t max_mag = 0;
    for (size_t i = 0; i <= constants::MAX_MAG; i++) {
        if (nfreqs->counts[i] != 0)
            max_mag = i;
    }
    /* first phase in scaling process, distribute out the
           last bucket, assume it is the smallest n(s) area, scale
           the rest by the same amount */
    auto bucket_size = uniq_vals_in_mag(max_mag, nfreqs->max_value);
    double C = 0.5 * bucket_size / nfreqs->counts[max_mag];
    for (size_t m = 0; m <= max_mag; m++) {
        bucket_size = uniq_vals_in_mag(m, nfreqs->max_value);
        nfreqs->counts[m] = 0.5 + nfreqs->counts[m] * C / bucket_size;
        if (table->counts[m] != 0 && nfreqs->counts[m] < 1) {
            nfreqs->counts[m] = 1;
        }
    }
    // print_mag_table(nfreqs, "first_phase");
    /* second step in scaling process, to make the first freq
           less than or equal to TOPFREQ
        */
    if (nfreqs->counts[0] > constants::TOPFREQ) {
        C = 1.0 * constants::TOPFREQ / nfreqs->counts[0];
        nfreqs->counts[0] = constants::TOPFREQ;
        /* scale all the others, rounding up so not zero anywhere,
               and at the same time, spread right across the bucketed
               range
            */
        for (uint8_t m = 1; m <= max_mag; m++) {
            nfreqs->counts[m] = 0.5 + nfreqs->counts[m] * C;
            if (table->counts[m] != 0 && nfreqs->counts[m] == 0) {
                nfreqs->counts[m] = 1;
            }
        }
    }
    // print_mag_table(nfreqs, "second_phase");

    /* now, what does it all add up to? */
    uint64_t M = 0;
    for (size_t m = 0; m <= max_mag; m++) {
        M += nfreqs->counts[m] * uniq_vals_in_mag(m, nfreqs->max_value);
    }
    /* fourth phase, round up to a power of two and then redistribute */
    uint64_t target_power = next_power_of_two(M);
    uint64_t excess = target_power - M;
    /* flow that excess count backwards to the beginning of
           the selectors array, spreading it out across the buckets...
        */
    for (int8_t m = int8_t(max_mag); m >= 0; m--) {
        double ratio = 1.0 * excess / M;
        uint64_t adder = ratio * nfreqs->counts[m];
        excess -= uniq_vals_in_mag(m, nfreqs->max_value) * adder;
        M -= uniq_vals_in_mag(m, nfreqs->max_value) * nfreqs->counts[m];
        nfreqs->counts[m] += adder;
    }
    if (excess != 0) {
        for (size_t m = 0; m <= max_mag; m++) {
            if (nfreqs->counts[m] != 0) {
                nfreqs->counts[m] += excess;
                break;
            }
        }
    }
    print_mag_table(nfreqs, "final_phase");

    M = 0;
    for (size_t i = 0; i <= max_mag; i++) {
        M += int64_t(nfreqs->counts[i] * uniq_vals_in_mag(i, nfreqs->max_value));
    }

    if (!is_power_of_two(M)) {
        fprintf(stderr, "ERROR! not power of 2 after normalization = %lu\n", M);
        exit(EXIT_FAILURE);
    }

    return nfreqs;
}
}