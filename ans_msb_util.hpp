#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "ans_util.hpp"

//#define ANS_DEBUG 1

namespace ans_msb {

namespace constants {
    constexpr uint32_t MAX_VAL = 1024;
    constexpr uint32_t M = MAX_VAL * 16;
    constexpr uint32_t MASK_M = M - 1;
    constexpr uint32_t LOG2_M = log2(M);
    constexpr uint8_t OUTPUT_BASE_LOG2 = 32;
    constexpr uint64_t OUTPUT_BASE = 1ULL << OUTPUT_BASE_LOG2;
    constexpr uint64_t NORM_LOWER_BOUND = 1ULL << 31;
}

using counts = uint64_t[constants::MAX_VAL + 1];

void print_counts(const counts& tb, std::string name)
{
    std::cout << "COUNTS = (";
    for (size_t i = 0; i <= constants::MAX_VAL; i++) {
        if (tb[i] != 0)
            std::cout << i << "=" << tb[i] << "\n";
    }
    std::cout << ")" << std::endl;
}

#pragma pack(1)
struct enc_table_entry {
    uint64_t freq;
    uint64_t base;
    uint64_t SUB;
};

std::ostream& operator<<(std::ostream& os, const enc_table_entry& t)
{
    os << "<f=" << t.freq << ",b=" << t.base << ",SUB=" << t.SUB << ">";
    return os;
}

using enc_model = enc_table_entry[constants::MAX_VAL + 1];

#pragma pack(1)
struct dec_table_entry {
    uint16_t freq;
    uint16_t except_bytes;
    uint64_t offset;
    uint32_t mapped_num;
    uint64_t except_mask;
};

using dec_model = dec_table_entry[constants::M];

uint16_t mapping(uint32_t x)
{
    uint8_t lzb = 3 - (__builtin_clz(x) >> 3);
    return (x >> (lzb << 3)) + (lzb << 8);
}

uint32_t undo_mapping(uint32_t x)
{
    if (x <= 256)
        return x;
    if (x <= 512)
        return ((x - 256) << 8);
    if (x <= 768)
        return ((x - 512) << 16);
    return ((x - 768) << 24);
}

uint32_t exception_bytes(uint32_t x)
{
    if (x <= 256)
        return 0;
    if (x <= 512)
        return 1;
    if (x <= 768)
        return 2;
    return 3;
}

uint32_t undo_mapping(const dec_table_entry& entry, const uint8_t*& except_ptr)
{
    auto except_u32 = reinterpret_cast<const uint32_t*>(except_ptr);
    uint64_t except_u64 = (*except_u32);
    uint32_t num = entry.mapped_num + (except_u64 & entry.except_mask);
    except_ptr += entry.except_bytes;
    return num;
}

uint16_t mapping_and_exceptions(uint32_t x, uint8_t*& except_out)
{
    if (x <= 256)
        return x;
    if (x <= (1 << 16)) {
        // one exception byte
        // std::cout << "x = " << x << " exception byte = " << (x & 0xFF) << std::endl;
        *except_out-- = x & 0xFF;
        return (x >> 8) + 256;
    }
    if (x <= (1 << 24)) {
        // two exception byte
        *except_out-- = (x >> 8) & 0xFF;
        *except_out-- = x & 0xFF;
        return (x >> 16) + 512;
    }
    // three exception byte
    *except_out-- = (x >> 16) & 0xFF;
    *except_out-- = (x >> 8) & 0xFF;
    *except_out-- = x & 0xFF;
    return (x >> 24) + 768;
}

uint16_t mapping_alistair(uint32_t x)
{
    if (x <= 256)
        return x;
    if (x <= (1 << 16))
        return (x >> 8) + 256;
    if (x <= (1 << 24))
        return (x >> 16) + 512;
    return (x >> 24) + 768;
}

std::vector<uint64_t> normalize_freqs(const counts& freqs, size_t target_power)
{
    // print_counts(freqs, "initital");
    std::vector<uint64_t> nfreqs(freqs, freqs + constants::MAX_VAL + 1);
    uint32_t n = 0;
    uint64_t initial_sum = 0;
    for (size_t i = 1; i < nfreqs.size(); i++) {
        if (freqs[i] != 0) {
            n = i + 1;
            initial_sum += freqs[i];
        }
    }
    /* first phase in scaling process, distribute out the
       last bucket, assume it is the smallest n(s) area, scale
       the rest by the same amount */
    double C = double(target_power) / double(initial_sum);
    for (size_t i = 1; i < n; i++) {
        nfreqs[i] = 0.95 * nfreqs[i] * C;
        if (freqs[i] != 0 && nfreqs[i] < 1) {
            nfreqs[i] = 1;
        }
    }

    /* now, what does it all add up to? */
    uint64_t M = 0;
    for (size_t m = 0; m < n; m++) {
        M += nfreqs[m];
    }
    /* fourth phase, round up to a power of two and then redistribute */
    uint64_t excess = target_power - M;
    /* flow that excess count backwards to the beginning of
       the selectors array, spreading it out across the buckets...
    */
    for (int64_t m = int64_t(n - 1); m >= 1; m--) {
        double ratio = double(excess) / double(M);
        uint64_t adder = ratio * nfreqs[m];
        if (adder > excess) {
            adder = excess;
        }
        excess -= adder;
        M -= nfreqs[m];
        nfreqs[m] += adder;
    }
    if (excess != 0) {
        nfreqs[0] += excess;
    }

    M = 0;
    for (size_t i = 0; i < n; i++) {
        M += nfreqs[i];
    }
    if (!ans::is_power_of_two(M)) {
        fprintf(stderr, "ERROR! not power of 2 after normalization = %lu", M);
        exit(EXIT_FAILURE);
    }
    // print_counts(freqs, "final");
    return nfreqs;
}

bool create_enc_model(std::vector<uint8_t>& enc_models, const counts& cnts)
{
    // (0) if all is 0 do nothing
    if (std::all_of(cnts, cnts + constants::MAX_VAL + 1,
            [](uint64_t i) { return i == 0; })) {
        return true;
    }

    // (1) normalize the counts
    auto norm_counts = ans_msb::normalize_freqs(cnts, constants::M);

    // (2) create the encoding model
    size_t model_size = sizeof(ans_msb::enc_model);
    std::vector<uint8_t> new_model(model_size);
    auto model_ptr = reinterpret_cast<ans_msb::enc_model*>(new_model.data());
    auto& model = *model_ptr;

    // (2a) fill the tables
    uint64_t base = 0;
    for (size_t i = 0; i <= constants::MAX_VAL; i++) {
        model[i].freq = norm_counts[i];
        model[i].base = base;
        base += norm_counts[i];
    }
    const uint64_t tmp = ((constants::NORM_LOWER_BOUND / constants::M) * ans_msb::constants::OUTPUT_BASE);
    for (size_t k = 0; k <= constants::MAX_VAL; k++) {
        model[k].SUB = tmp * model[k].freq;
    }
    enc_models.insert(enc_models.end(), new_model.begin(), new_model.end());
    return false;
}

void create_dec_model(std::vector<uint8_t>& dec_models, const ans_msb::enc_model& enc_model)
{
    // (1) determine model size
    size_t model_size = sizeof(ans_msb::dec_model);
    std::vector<uint8_t> new_model(model_size);
    auto model_ptr = reinterpret_cast<ans_msb::dec_model*>(new_model.data());
    auto& model = *model_ptr;

    // (2) create csum table for decoding
    size_t base = 0;
    for (size_t j = 0; j <= constants::MAX_VAL; j++) {
        auto cur_freq = enc_model[j].freq;
        for (size_t k = 0; k < cur_freq; k++) {
            model[base + k].mapped_num = undo_mapping(j);
            model[base + k].except_bytes = exception_bytes(j);
            model[base + k].except_mask = (1ULL << (model[base + k].except_bytes << 3)) - 1;
            model[base + k].freq = cur_freq;
            model[base + k].offset = k;
        }
        base += cur_freq;
    }
    dec_models.insert(dec_models.end(), new_model.begin(), new_model.end());
}

uint64_t encode_num(const enc_model& model, uint64_t state, uint32_t num, uint8_t*& out)
{
    const auto& entry = model[num];
    // (1) normalize
    if (state >= entry.SUB) {
        ans::output_unit<ans::constants::OUTPUT_BASE_LOG2>(out, state);
    }
    // (2) transform state
    uint64_t next = ((state / entry.freq) * ans_msb::constants::M) + (state % entry.freq) + entry.base;
    return next;
}

const dec_table_entry& decode_num(const dec_model& model, uint64_t& state, const uint8_t*& in, size_t& enc_size)
{
    uint64_t state_mod_M = state & ans_msb::constants::MASK_M;
    const auto& entry = model[state_mod_M];
    state = entry.freq * (state >> ans_msb::constants::LOG2_M) + entry.offset;
    if (enc_size && state < ans_msb::constants::NORM_LOWER_BOUND) {
        ans::input_unit<ans::constants::OUTPUT_BASE_LOG2>(in, state, enc_size);
    }
    return entry;
}
}