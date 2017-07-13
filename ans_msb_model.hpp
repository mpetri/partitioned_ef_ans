#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "ans_msb_util.hpp"
#include "ans_util.hpp"

struct block_header {
    uint32_t model_id = 0;
    uint32_t final_state_bytes = 0;
    uint32_t num_ans_u32s = 0;
};

struct model_max_1d {
    static const uint32_t NUM_MODELS = 16;
    static const uint8_t MAX_MAG = 32;
    static constexpr std::array<uint32_t, NUM_MODELS> SEL2MAG{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12,
        14, 16, 19, 22, 32 };
    static constexpr std::array<uint32_t, MAX_MAG + 1> MAG2SEL{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
        10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
    static uint32_t pick_model(uint32_t const* in, size_t n)
    {
        uint32_t max_val = 0;
        for (size_t i = 0; i < n; i++) {
            max_val = std::max(max_val, in[i] + 1);
        }
        uint32_t max_mag = ans::magnitude(max_val);
        return MAG2SEL[max_mag];
    }

    static void write_block_header(const block_header& bh, std::vector<uint8_t>& out)
    {
    }

    static void read_block_header(block_header& bh, uint8_t const*& in)
    {
    }
};

namespace quasi_succinct {

template <typename model_type>
struct ans_msb_model {
    static const uint64_t block_size = ans::constants::BLOCK_SIZE;
    static const uint32_t NUM_MODELS = model_type::NUM_MODELS;
    using ans_msb_counts = ans_msb::counts[NUM_MODELS];

    static std::vector<uint8_t> create_empty_counts()
    {
        size_t count_size = sizeof(ans_msb::counts);
        return std::vector<uint8_t>(count_size, 0);
    }

    static std::vector<uint8_t> create_enc_model_from_counts(const std::vector<uint8_t>& cntsu8)
    {
        auto counts_ptr = reinterpret_cast<const ans_msb_counts*>(cntsu8.data());
        const ans_msb_counts& counts = *counts_ptr;
        size_t pointers_to_models = NUM_MODELS * sizeof(uint64_t);
        std::vector<uint8_t> enc_models(pointers_to_models, 0);
        for (size_t i = 0; i < NUM_MODELS; i++) {
            size_t model_offset = enc_models.size();
            // (1) create the model
            bool empty_model = ans_msb::create_enc_model(enc_models, counts[i]);
            auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(enc_models.data()) + i;
            // (1) store offset of model data in byte stream
            if (empty_model) {
                *model_offset_u64_ptr = 0;
            } else {
                *model_offset_u64_ptr = model_offset;
            }
        }
        return enc_models;
    }

    static std::vector<uint8_t> create_dec_models(const std::vector<uint8_t>& enc_models_u8)
    {
        auto enc_models = reinterpret_cast<const uint64_t*>(enc_models_u8.data());
        size_t pointers_to_models = NUM_MODELS * sizeof(uint64_t);
        std::vector<uint8_t> dec_models_u8(pointers_to_models, 0);
        for (size_t i = 0; i < NUM_MODELS; i++) {
            size_t dec_model_offset = dec_models_u8.size();
            if (enc_models[i] != 0) {
                size_t enc_model_offset = enc_models[i];
                auto enc_model_ptr = reinterpret_cast<const ans_msb::enc_model*>(enc_models_u8.data() + enc_model_offset);
                size_t model_size = dec_models_u8.size();
                ans_msb::create_dec_model(dec_models_u8, *enc_model);
                model_size = dec_models_u8.size() - model_size;
                std::cout << "CREATE_DEC_MODEL(" << i << ") bytes = " << model_size << std::endl;
                auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(dec_models_u8.data()) + i;
                *model_offset_u64_ptr = dec_model_offset;
            } else {
                auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(dec_models_u8.data()) + i;
                *model_offset_u64_ptr = 0;
            }
        }
        return dec_models_u8;
    }

    static void
    model(std::vector<uint8_t>& cntsu8, uint32_t const* in, uint32_t /*sum_of_values*/, size_t n)
    {
        auto counts_ptr = reinterpret_cast<ans_msb_counts*>(cntsu8.data());
        auto& counts = *counts_ptr;
        auto model_id = model_type::pick_model(in, n);
        for (size_t i = 0; i < n; i++) {
            auto msb_val = ans_msb::mapping(in[i] + 1);
            counts[model_id][msb_val]++;
        }
    }

    static void encode(uint32_t const* in, uint32_t /*sum_of_values*/,
        size_t n, std::vector<uint8_t>& out, const std::vector<uint8_t>& enc_model_u8)
    {
    }

    static uint8_t const*
    decode(uint8_t const* in, uint32_t* out,
        uint32_t /* sum_of_values */, size_t n, uint8_t const* dec_model_u8)
    {

        return in;
    }
}
};
}
