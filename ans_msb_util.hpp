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
    constexpr uint32_t FRAME_SIZE_FACTOR = 16;
    constexpr uint32_t MAX_VAL = 1024;
    constexpr uint32_t MAX_M = MAX_VAL * FRAME_SIZE_FACTOR;
    constexpr uint8_t OUTPUT_BASE_LOG2 = 32;
    constexpr uint64_t OUTPUT_BASE = 1ULL << OUTPUT_BASE_LOG2;
    constexpr uint64_t NORM_LOWER_BOUND = 1ULL << 31;
    constexpr uint64_t VBYTE_THRESHOLD = 10;
    constexpr uint64_t COMPACT_THRESHOLD = 16;
}

using counts = uint64_t[constants::MAX_VAL + 1];

template <class t_vec>
void print_cnts(const t_vec& cnts)
{
    for (size_t i = 0; i <= constants::MAX_VAL; i++) {
        if (cnts[i] != 0) {
            std::cout << i << " = " << cnts[i] << std::endl;
        }
    }
}

std::pair<double, uint64_t>
compute_entropy(const counts& cnts)
{
    double H = 0;
    double N = 0;
    for (size_t i = 0; i <= constants::MAX_VAL; i++)
        N += cnts[i];
    if (N == 0)
        return { 0, 0 };
    for (size_t i = 0; i <= constants::MAX_VAL; i++) {
        double n = cnts[i];
        if (cnts[i])
            H += n * log2(N / n);
    }
    return { H, N };
}

double
compute_loss(counts& A, counts& B, std::pair<double, uint64_t> H_a, std::pair<double, uint64_t> H_b)
{
    if (H_a.second == 0 || H_b.second == 0) // cant merge empty stuff
        return std::numeric_limits<double>::max();

    // compute combined entropy
    double Nc = H_a.second + H_b.second;
    double Hc = 0.0;
    for (size_t i = 0; i <= constants::MAX_VAL; i++) {
        if (A[i] || B[i]) {
            Hc += double(A[i] + B[i]) * log2(Nc / double(A[i] + B[i]));
        }
    }
    return Hc - H_a.first - H_b.first;
}

template <class t_cnts>
std::pair<uint64_t, uint64_t>
find_min_pair(t_cnts& cnts, size_t num_models, std::vector<std::pair<double, uint64_t>>& model_entropy)
{
    double min_loss = std::numeric_limits<double>::max();
    std::pair<uint64_t, uint64_t> min_pair{ 0, 0 };
    for (size_t i = 0; i < num_models; i++) {
        for (size_t j = i + 1; j < num_models; j++) {
            double pair_loss = compute_loss(cnts[i], cnts[j], model_entropy[i], model_entropy[j]);
            if (pair_loss < min_loss) {
                min_loss = pair_loss;
                min_pair.first = i;
                min_pair.second = j;
            }
        }
    }
    std::cout << "find_min_pair = <" << min_pair.first << "," << min_pair.second << "> loss = " << min_loss << std::endl;
    return min_pair;
}

template <class t_cnts>
void merge_models(t_cnts& counts, std::vector<std::pair<double, uint64_t>>& model_entropy, uint64_t from, uint64_t to)
{
    model_entropy[from] = { 0.0, 0 };
    for (size_t i = 0; i <= constants::MAX_VAL; i++) {
        counts[to][i] += counts[from][i];
        counts[from][i] = 0;
    }
    std::cout << "merge_moddels(from=" << from << ",to=" << to << ")" << std::endl;
    model_entropy[to] = compute_entropy(counts[to]);
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

struct enc_model {
    uint64_t M;
    enc_table_entry table[0];
};

#pragma pack(1)
struct dec_table_entry {
    uint16_t freq;
    uint16_t except_bytes;
    uint64_t offset;
    uint32_t mapped_num;
    uint64_t except_mask;
};

struct dec_model {
    uint64_t M;
    uint64_t MASK_M;
    uint64_t LOG2_M;
    dec_table_entry table[0];
};

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
    uint64_t M = std::numeric_limits<uint64_t>::max();
    float fudge = 1.0;
    while (M > target_power) {
        fudge -= 0.01;
        for (size_t i = 1; i < n; i++) {
            nfreqs[i] = fudge * freqs[i] * C;
            if (freqs[i] != 0 && nfreqs[i] < 1) {
                nfreqs[i] = 1;
            }
        }

        /* now, what does it all add up to? */
        M = 0;
        for (size_t m = 0; m < n; m++) {
            M += nfreqs[m];
        }
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
        std::cerr << "ERROR! not power of 2 after normalization = " << M << std::endl;
        exit(EXIT_FAILURE);
    }
    return nfreqs;
}

bool create_enc_model(std::vector<uint8_t>& enc_models, const counts& cnts)
{
    // (0) if all is 0 do nothing
    if (std::all_of(cnts, cnts + constants::MAX_VAL + 1,
            [](uint64_t i) { return i == 0; })) {
        return true;
    }

    // (1) find a target power of 2
    uint64_t uniq_syms = 0;
    for (size_t i = 0; i < constants::MAX_VAL + 1; i++) {
        if (cnts[i] != 0)
            uniq_syms++;
    }
    uint64_t target_M = uniq_syms * constants::FRAME_SIZE_FACTOR;
    uint64_t target_PTWO = target_M;
    if (!ans::is_power_of_two(target_PTWO))
        target_PTWO = ans::next_power_of_two(target_PTWO);

    // (1) normalize the counts
    // std::cout << "unique syms = " << uniq_syms << std::endl;
    // std::cout << "actual counts = " << std::endl;
    // print_cnts(cnts);
    auto norm_counts = ans_msb::normalize_freqs(cnts, target_PTWO);
    // std::cout << "normalized counts = " << std::endl;
    // print_cnts(norm_counts);

    // (2) create the encoding model
    size_t model_size = sizeof(ans_msb::enc_model) + (constants::MAX_VAL + 1) * sizeof(enc_table_entry);
    std::vector<uint8_t> new_model(model_size);
    auto model_ptr = reinterpret_cast<ans_msb::enc_model*>(new_model.data());
    auto& model = *model_ptr;
    model.M = target_PTWO;
    std::cout << "M = " << model.M << std::endl;

    // (2a) fill the tables
    uint64_t base = 0;
    for (size_t i = 0; i <= constants::MAX_VAL; i++) {
        model.table[i].freq = norm_counts[i];
        model.table[i].base = base;
        base += norm_counts[i];
    }
    const uint64_t tmp = ((constants::NORM_LOWER_BOUND / model.M) * ans_msb::constants::OUTPUT_BASE);
    for (size_t k = 0; k <= constants::MAX_VAL; k++) {
        model.table[k].SUB = tmp * model.table[k].freq;
    }
    enc_models.insert(enc_models.end(), new_model.begin(), new_model.end());
    return false;
}

void create_dec_model(std::vector<uint8_t>& dec_models, const ans_msb::enc_model& enc_model)
{
    // (1) determine model size
    size_t model_size = sizeof(ans_msb::dec_model) + enc_model.M * sizeof(dec_table_entry);
    std::vector<uint8_t> new_model(model_size);
    auto model_ptr = reinterpret_cast<ans_msb::dec_model*>(new_model.data());
    auto& model = *model_ptr;

    model.M = enc_model.M;
    model.LOG2_M = log2(enc_model.M);
    model.MASK_M = model.M - 1;

    // (2) create csum table for decoding
    size_t base = 0;
    for (size_t j = 0; j <= constants::MAX_VAL; j++) {
        auto cur_freq = enc_model.table[j].freq;
        for (size_t k = 0; k < cur_freq; k++) {
            model.table[base + k].mapped_num = undo_mapping(j);
            model.table[base + k].except_bytes = exception_bytes(j);
            model.table[base + k].except_mask = (1ULL << (model.table[base + k].except_bytes << 3)) - 1;
            model.table[base + k].freq = cur_freq;
            model.table[base + k].offset = k;
        }
        base += cur_freq;
    }
    dec_models.insert(dec_models.end(), new_model.begin(), new_model.end());
}

uint64_t encode_num(const enc_model& model, uint64_t state, uint32_t num, uint8_t*& out)
{
    const auto& entry = model.table[num];
    // (1) normalize
    if (state >= entry.SUB) {
        ans::output_unit<ans::constants::OUTPUT_BASE_LOG2>(out, state);
    }
    // (2) transform state
    uint64_t next = ((state / entry.freq) * model.M) + (state % entry.freq) + entry.base;
    return next;
}

const dec_table_entry& decode_num(const dec_model& model, uint64_t& state, const uint8_t*& in, size_t& enc_size)
{
    uint64_t state_mod_M = state & model.MASK_M;
    const auto& entry = model.table[state_mod_M];
    state = entry.freq * (state >> model.LOG2_M) + entry.offset;
    if (enc_size && state < ans_msb::constants::NORM_LOWER_BOUND) {
        ans::input_unit<ans::constants::OUTPUT_BASE_LOG2>(in, state, enc_size);
    }
    return entry;
}
}
