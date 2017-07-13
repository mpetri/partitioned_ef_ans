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
        // (1) determine and encode model id
        block_header bh;
        bh.model_id = model_type::pick_model(in, n);

        if (bh.model_id == 0) { // all 1s. continue
            model_type::write_block_header(bh, out);
            return;
        }

        static std::array<uint8_t, block_size * 8> tmp_out_buf;
        static std::array<uint8_t, block_size * 8> exception_out_buf;

        // (2) reverse encode the block using the selected ANS model
        auto model_ptrs = reinterpret_cast<const uint64_t*>(enc_model_u8.data());
        size_t model_offset = model_ptrs[bh.model_id];
        uint64_t state = 0;
        auto out_ptr = tmp_out_buf.data() + tmp_out_buf.size() - 1;
        auto except_ptr = exception_out_buf.data();
        auto out_start = out_ptr;
        auto except_start = except_ptr;
        auto cur_model = reinterpret_cast<const ans_msb::enc_model*>(enc_model_u8.data() + model_offset);
        for (size_t k = 0; k < n; k++) {
            uint32_t num = in[n - k - 1];
            uint32_t mapped_num = ans_msb::mapping_and_exceptions(num, except_ptr);
            state = ans_msb::encode_num(cur_model, state, mapped_num, out_ptr);
        }
        size_t enc_size = out_start - out_ptr;
        size_t u32s_written = enc_size / sizeof(uint32_t);

        // (3) write block header
        bh.final_state_bytes = ans::state_bytes(state);
        bh.num_ans_u32s = u32s_written;
        model_type::write_block_header(bh, out);

        // (4) write the final state
        ans::flush_state(state, out_ptr, bh.final_state_bytes);

        // (5) copy the ans output to the buffer
        size_t final_ans_size = out_start - out_ptr;
        out.insert(out.end(), out_ptr, out_ptr + final_ans_size);
        size_t final_except_size = except_ptr - except_start;
        out.insert(out.end(), except_start, except_start + final_except_size);
    }

    static uint8_t const*
    decode(uint8_t const* in, uint32_t* out,
        uint32_t /* sum_of_values */, size_t n, uint8_t const* dec_model_u8)
    {
        block_header bh;
        model_type::read_block_header(bh, in);

        // uniform block
        if (bh.model_id == 0) {
            for (size_t i = 0; i < n; i++)
                out[i] = 0;
            return in;
        }

        size_t enc_size = bh.num_ans_u32s * sizeof(uint32_t);
        auto model_ptrs = reinterpret_cast<const uint64_t*>(dec_model_u8);
        size_t model_offset = model_ptrs[bh.model_id];

        // (0) init the decoder
        uint64_t state = ans::init_decoder(in, bh.final_state_bytes);

        // (1) decode the ans parts
        auto cur_model = reinterpret_cast<const ans_msg::dec_model*>(dec_model_u8 + model_offset);
        bool has_exceptions = false;
        auto out_start = out;
        for (size_t k = 0; k < n; k++) {
            uint32_t dec_num = decode_num(cur_model, state, in, enc_size);
            if (dec_num > 256)
                has_exceptions = true;
            *out++ = dec_num; // substract one as OT has 0s and our compactest num is 1
        }

        // (2) deal with exceptions
        if (has_exceptions) {
            for (size_t k = 0; k < n; k++) {
                if (*out_start <= 256) {
                    out_start++;
                    continue;
                }
                if (*out_start > (1U << 8)) {
                    *out_start = (*out_start << 8) + *in++;
                }
                if (*out_start > (1U << 16)) {
                    *out_start = (*out_start << 8) + *in++;
                }
                if (*out_start > (1U << 24)) {
                    *out_start = (*out_start << 8) + *in++;
                }
                out_start++;
                continue;
            }
        }
        return in;
    }
}
};
}
