#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "ans_msb_util.hpp"
#include "ans_util.hpp"

struct msb_block_header {
    uint32_t model_id = 0;
    uint32_t final_state_bytes = 0;
    uint32_t num_ans_u32s = 0;
};

struct msb_model_max_1d {
    static const uint32_t NUM_MODELS = 16;
    static uint32_t pick_model(uint32_t const* in, size_t n)
    {
        static const std::vector<uint32_t> MAG2SEL{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
            10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
        uint32_t max_val = 0;
        for (size_t i = 0; i < n; i++) {
            max_val = std::max(max_val, in[i] + 1);
        }
        uint32_t max_mag = ans::magnitude(max_val);
        return MAG2SEL[max_mag];
    }

    static void write_block_header(const msb_block_header& bh, std::vector<uint8_t>& out)
    {
        uint8_t packed = ans::pack_two_4bit_nums(bh.model_id, bh.final_state_bytes);
        out.push_back(packed);
        if (bh.model_id != 0) {
            out.push_back(uint8_t(bh.num_ans_u32s));
        }
    }

    static void read_block_header(msb_block_header& bh, uint8_t const*& in)
    {
        uint8_t packed = *in++;
        auto model_and_fsb = ans::unpack_two_4bit_nums(packed);
        bh.model_id = model_and_fsb.first;
        bh.final_state_bytes = model_and_fsb.second;
        if (bh.model_id != 0) {
            bh.num_ans_u32s = *in++;
        }
    }
};

struct msb_model_minmax_2d {
    static const uint32_t NUM_MODELS = 16 * 16;
    static uint32_t pick_model(uint32_t const* in, size_t n)
    {
        static const std::vector<uint32_t> MAG2SEL{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
            10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
        uint32_t min_val = std::numeric_limits<uint32_t>::max();
        uint32_t max_val = 0;
        for (size_t i = 0; i < n; i++) {
            min_val = std::min(min_val, in[i] + 1);
            max_val = std::max(max_val, in[i] + 1);
        }
        uint32_t min_mag = ans::magnitude(min_val);
        uint32_t max_mag = ans::magnitude(max_val);
        return (MAG2SEL[max_mag] << 4) + MAG2SEL[min_mag];
    }

    static void write_block_header(const msb_block_header& bh, std::vector<uint8_t>& out)
    {
        out.push_back(bh.model_id);
        if (bh.model_id != 0) {
            out.push_back(bh.final_state_bytes);
            out.push_back(uint8_t(bh.num_ans_u32s));
        }
    }

    static void read_block_header(msb_block_header& bh, uint8_t const*& in)
    {
        bh.model_id = *in++;
        if (bh.model_id != 0) {
            bh.final_state_bytes = *in++;
            bh.num_ans_u32s = *in++;
        }
    }
};

struct msb_model_med90p_2d {
    static const uint32_t NUM_MODELS = 16 * 16;
    static uint32_t pick_model(uint32_t const* in, size_t n)
    {
        static const std::vector<uint32_t> MAG2SEL{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
            10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
        static std::vector<uint32_t> buf(ans::constants::BLOCK_SIZE);
        std::copy(in, in + n, buf.begin());
        std::sort(buf.begin(), buf.begin() + n);
        uint32_t mag_med = ans::magnitude(buf[n / 2] + 1);
        uint32_t mag_90p = ans::magnitude(buf[n * 0.9] + 1);
        return (MAG2SEL[mag_90p] << 4) + MAG2SEL[mag_med];
    }

    static void write_block_header(const msb_block_header& bh, std::vector<uint8_t>& out)
    {
        out.push_back(bh.model_id);
        if (bh.model_id != 0) {
            out.push_back(bh.final_state_bytes);
            out.push_back(uint8_t(bh.num_ans_u32s));
        }
    }

    static void read_block_header(msb_block_header& bh, uint8_t const*& in)
    {
        bh.model_id = *in++;
        if (bh.model_id != 0) {
            bh.final_state_bytes = *in++;
            bh.num_ans_u32s = *in++;
        }
    }
};

namespace quasi_succinct {

template <typename model_type>
struct ans_msb_model {
    static const uint64_t block_size = ans::constants::BLOCK_SIZE;
    static const uint32_t NUM_MODELS = model_type::NUM_MODELS;
    using ans_msb_counts_table = ans_msb::counts[NUM_MODELS];

    static std::vector<uint8_t> create_empty_counts()
    {
        size_t count_size = sizeof(ans_msb_counts_table);
        return std::vector<uint8_t>(count_size, 0);
    }

    static std::vector<uint8_t> create_enc_model_from_counts(const std::vector<uint8_t>& cntsu8)
    {
        auto counts_ptr = reinterpret_cast<const ans_msb_counts_table*>(cntsu8.data());
        auto& counts = *counts_ptr;
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
                ans_msb::create_dec_model(dec_models_u8, *enc_model_ptr);
                model_size = dec_models_u8.size() - model_size;
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
        auto counts_ptr = reinterpret_cast<ans_msb_counts_table*>(cntsu8.data());
        auto& counts = *counts_ptr;
        auto model_id = model_type::pick_model(in, n);
        for (size_t i = 0; i < n; i++) {
            auto msb_val = ans_msb::mapping_alistair(in[i] + 1);
            counts[model_id][msb_val]++;
        }
    }

    static void encode(uint32_t const* in, uint32_t /*sum_of_values*/,
        size_t n, std::vector<uint8_t>& out, const std::vector<uint8_t>& enc_model_u8)
    {
        // (1) determine and encode model id
        msb_block_header bh;
        bh.model_id = model_type::pick_model(in, n);

        if (bh.model_id == 0) { // all 1s. continue
            model_type::write_block_header(bh, out);
            return;
        }

        static std::vector<uint8_t> tmp_out_buf(block_size * 8);
        static std::vector<uint8_t> exception_out_buf(block_size * 8);

        // (2) reverse encode the block using the selected ANS model
        auto model_ptrs = reinterpret_cast<const uint64_t*>(enc_model_u8.data());
        size_t model_offset = model_ptrs[bh.model_id];
        uint64_t state = 0;
        auto out_ptr = tmp_out_buf.data() + tmp_out_buf.size() - 1;
        auto except_ptr = exception_out_buf.data() + exception_out_buf.size() - 1;
        auto out_start = out_ptr;
        auto except_start = except_ptr;
        auto cur_model = reinterpret_cast<const ans_msb::enc_model*>(enc_model_u8.data() + model_offset);
        for (size_t k = 0; k < n; k++) {
            uint32_t num = in[n - k - 1] + 1;
            uint32_t mapped_num = ans_msb::mapping_and_exceptions(num, except_ptr);
            state = ans_msb::encode_num(*cur_model, state, mapped_num, out_ptr);
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
        size_t final_except_size = except_start - except_ptr;
        if (final_except_size)
            out.insert(out.end(), except_ptr + 1, except_ptr + final_except_size + 1);
    }

    static uint8_t const*
    decode(uint8_t const* in, uint32_t* out,
        uint32_t /* sum_of_values */, size_t n, uint8_t const* dec_model_u8)
    {
        msb_block_header bh;
        model_type::read_block_header(bh, in);

        // uniform block
        if (bh.model_id == 0) {
            for (size_t i = 0; i < n; i++)
                out[i] = 0;
            return in;
        }

        size_t ans_enc_size = bh.num_ans_u32s * sizeof(uint32_t);
        auto model_ptrs = reinterpret_cast<const uint64_t*>(dec_model_u8);
        size_t model_offset = model_ptrs[bh.model_id];

        // (0) init the decoder
        uint64_t state = ans::init_decoder(in, bh.final_state_bytes);

        // (1) decode the ans parts and the exceptions at the same time
        auto except_ptr = in + ans_enc_size;
        auto cur_model = reinterpret_cast<const ans_msb::dec_model*>(dec_model_u8 + model_offset);
        for (size_t k = 0; k < n; k++) {
            const auto& dec_entry = decode_num(*cur_model, state, in, ans_enc_size);
            *out++ = ans_msb::undo_mapping(dec_entry, except_ptr) - 1;
        }
        return except_ptr;
    }
};
}
