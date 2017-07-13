#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

//#define ANS_DEBUG 1

namespace ans_msb {

namespace constants {
    const uint32_t MAX_VAL = 1024;
    const uint32_t FRAME_SIZE = 1024 * 16;
    const uint8_t OUTPUT_BASE_LOG2 = 32;
    const uint64_t OUTPUT_BASE = 1ULL << OUTPUT_BASE_LOG2;
    const uint64_t NORM_LOWER_BOUND = 1ULL << 24;
}

using counts = uint64_t[constants::MAX_VAL + 1];

#pragma pack(1)
struct enc_table_entry {
    uint64_t freq;
    uint64_t base;
    uint64_t SUB;
};

using enc_model = enc_table_entry[constants::MAX_VAL];

#pragma pack(1)
struct dec_table_entry {
    uint64_t freq;
    uint64_t base;
    uint64_t SUB;
};

using dec_model = dec_table_entry[constants::FRAME_SIZE];

uint16_t mapping(uint32_t x)
{
    uint8_t lzb = 3 - (__builtin_clz(x - 1) >> 3);
    return (x >> (lzb << 3)) + (lzb << 8);
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
    if (!is_power_of_two(M)) {
        quit("ERROR! not power of 2 after normalization = %lu", M);
    }

    return nfreqs;
}

bool create_enc_model(std::vector<uint8_t>& enc_models, counts& cnts)
{
    // (0) if all is 0 do nothing
    if (std::all_of(cnts, cnts + constants::MAX_VAL + 1,
            [](uint64_t i) { return i == 0; })) {
        return true;
    }

    // (1) normalize the counts
    auto norm_counts = ans_msb::normalize_freqs(cnts, constants::FRAME_SIZE);

    // (2) create the encoding model
    size_t model_size = sizeof(ans_msb::enc_model);
    std::vector<uint8_t> new_model(model_size);
    auto model_ptr = reinterpret_cast<ans_msb::enc_model*>(new_model.data());
    auto& model = *model_ptr;

    // (2a) fill the tables
    uint64_t cumsum = 0;
    uint64_t j = 0;
    for (size_t i = 0; i <= constants::MAX_VAL; i++) {
        if (norm_counts[i] == 0)
            continue;
        for (size_t j = 0; j < norm_counts[i]; j++) {
            model[j].freq = norm_counts[i];
            model[j].base = cumsum;
        }
        cumsum += norm_counts[i];
    }
    for (size_t k = 0; k <= constants::MAX_VAL; k++) {
        model[k].SUB = ((constants::NORM_LOWER_BOUND / constants::FRAME_SIZE) * ans_msb::constants::OUTPUT_BASE)
            * model[k].freq;
    }
    enc_models.insert(enc_models.end(), new_model.begin(), new_model.end());
    return false;
}

static void create_dec_model(std::vector<uint8_t>& dec_models, const ans_msb::enc_model& enc_model)
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
            model[base + k].sym = j;
            model[base + k].freq = cur_freq;
            model[base + k].offset = k;
        }
        base += cur_freq;
    }
    dec_models.insert(dec_models.end(), new_model.begin(), new_model.end());
}
}