#pragma once

#include "ans_packed_util.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

struct block_header {
    uint32_t model_id = 0;
    uint32_t final_state_bytes = 0;
    uint32_t num_ans_u32s = 0;
};

struct model_max_1d {
    static const uint32_t NUM_MODELS = 16;
    static uint8_t pick_model(uint32_t const* in, size_t n)
    {
        uint32_t max_val = 0;
        for (size_t i = 0; i < n; i++) {
            max_val = std::max(max_val, in[i] + 1);
        }
        uint8_t max_mag = ans_packed::magnitude(max_val);
        return ans_packed::constants::MAG2SEL[max_mag];
    }

    static void write_block_header(const block_header& bh, std::vector<uint8_t>& out)
    {
        uint8_t packed = ans_packed::pack_two_4bit_nums(bh.model_id, bh.final_state_bytes);
        out.push_back(packed);
        if (bh.model_id != 0) {
            out.push_back(uint8_t(bh.num_ans_u32s));
        }
    }

    static void read_block_header(block_header& bh, uint8_t const*& in)
    {
        uint8_t packed = *in++;
        auto model_and_fsb = ans_packed::unpack_two_4bit_nums(packed);
        bh.model_id = model_and_fsb.first;
        bh.final_state_bytes = model_and_fsb.second;
        if (bh.model_id != 0) {
            bh.num_ans_u32s = *in++;
        }
    }
};

struct model_minmax_2d {
    static const uint32_t NUM_MODELS = 16 * 16;
    static uint8_t pick_model(uint32_t const* in, size_t n)
    {
        uint32_t max_val = 0;
        uint32_t min_val = std::numeric_limits<uint32_t>::max();
        for (size_t i = 0; i < n; i++) {
            min_val = std::min(min_val, in[i] + 1);
            max_val = std::max(max_val, in[i] + 1);
        }
        uint8_t min_mag = ans_packed::magnitude(min_val);
        uint8_t max_mag = ans_packed::magnitude(max_val);
        return (ans_packed::constants::MAG2SEL[max_mag] << 4) + ans_packed::constants::MAG2SEL[min_mag];
    }

    static void write_block_header(const block_header& bh, std::vector<uint8_t>& out)
    {
        out.push_back(uint8_t(bh.model_id));
        if (bh.model_id != 0) {
            out.push_back(uint8_t(bh.final_state_bytes));
            out.push_back(uint8_t(bh.num_ans_u32s));
        }
    }

    static void read_block_header(block_header& bh, uint8_t const*& in)
    {
        bh.model_id = *in++;
        if (bh.model_id != 0) {
            bh.final_state_bytes = *in++;
            bh.num_ans_u32s = *in++;
        }
    }
};

struct model_med90p_2d {
    static const uint32_t NUM_MODELS = 16 * 16;
    static uint8_t pick_model(uint32_t const* in, size_t n)
    {
        static std::vector<uint32_t> tmp(ans_packed::constants::BLOCK_SIZE);
        std::copy(in, in + n, tmp.begin());
        std::sort(tmp.begin(), tmp.begin() + n);
        uint32_t val_med = tmp[n / 2] + 1;
        uint32_t val_90p = tmp[0.9 * n] + 1;
        uint8_t mag_med = ans_packed::magnitude(val_med);
        uint8_t mag_90p = ans_packed::magnitude(val_90p);
        if (tmp[n - 1] + 1 == 1) {
            return 0;
        }
        uint8_t model_id = (ans_packed::constants::MAG2SEL[mag_90p] << 4) + ans_packed::constants::MAG2SEL[mag_med];
        if (model_id == 0)
            return 1;
        return model_id;
    }

    static void write_block_header(const block_header& bh, std::vector<uint8_t>& out)
    {
        out.push_back(uint8_t(bh.model_id));
        if (bh.model_id != 0) {
            out.push_back(uint8_t(bh.final_state_bytes));
            out.push_back(uint8_t(bh.num_ans_u32s));
        }
    }

    static void read_block_header(block_header& bh, uint8_t const*& in)
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
struct ans_packed_model {
    static const uint64_t block_size = ans_packed::constants::BLOCK_SIZE;
    static const uint32_t NUM_MODELS = model_type::NUM_MODELS;
    using ans_packed_counts = ans_packed::mag_table[NUM_MODELS];

    static std::vector<uint8_t> create_empty_counts()
    {
        size_t count_size = sizeof(ans_packed_counts);
        return std::vector<uint8_t>(count_size, 0);
    }

    static size_t set_compact_model(size_t model_offset)
    {
        return model_offset | 0x8000000000000000ULL;
    }

    static bool is_compact_model(size_t& model_offset)
    {
        if (model_offset & 0x8000000000000000ULL) {
            model_offset ^= 0x8000000000000000ULL;
            return true;
        }
        return false;
    }

    static std::vector<uint8_t> create_enc_model_from_counts(const std::vector<uint8_t>& cntsu8)
    {
        auto counts_ptr = reinterpret_cast<const ans_packed_counts*>(cntsu8.data());
        const ans_packed_counts& counts = *counts_ptr;
        size_t pointers_to_models = NUM_MODELS * sizeof(uint64_t);
        std::vector<uint8_t> enc_models(pointers_to_models, 0);

        for (size_t i = 0; i < NUM_MODELS; i++) {
            // (1) store offset of model data in byte stream
            size_t model_offset = enc_models.size();
            // (2) create the model

            bool empty_model = false;
            bool compact_model = false;
            if (counts[i].max_value < ans_packed::constants::COMPACT_MODEL_THRESHOLD) {
                empty_model = create_enc_model(enc_models, counts[i]);
            } else {
                compact_model = true;
                empty_model = create_enc_model_compact(enc_models, counts[i]);
                std::cout << "CREATE_COMPACT_MODEL(" << i << ") max_value = " << counts[i].max_value << std::endl;
            }
            auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(enc_models.data()) + i;
            if (empty_model) {
                *model_offset_u64_ptr = 0;
            } else {
                if (compact_model)
                    *model_offset_u64_ptr = set_compact_model(model_offset);
                else
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
                if (is_compact_model(enc_model_offset)) {
                    auto enc_model_ptr = reinterpret_cast<const ans_packed::enc_model_compact*>(enc_models_u8.data() + enc_model_offset);
                    const ans_packed::enc_model_compact& enc_model = *enc_model_ptr;
                    create_dec_model_compact(dec_models_u8, enc_model);
                    std::cout << "CREATE_COMPACT_DEC_MODEL(" << i << ")" << std::endl;
                    auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(dec_models_u8.data()) + i;
                    *model_offset_u64_ptr = set_compact_model(dec_model_offset);
                } else {
                    auto enc_model_ptr = reinterpret_cast<const ans_packed::enc_model*>(enc_models_u8.data() + enc_model_offset);
                    const ans_packed::enc_model& enc_model = *enc_model_ptr;
                    create_dec_model(dec_models_u8, enc_model);
                    auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(dec_models_u8.data()) + i;
                    *model_offset_u64_ptr = dec_model_offset;
                }

            } else {
                auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(dec_models_u8.data()) + i;
                *model_offset_u64_ptr = 0;
            }
        }
        return dec_models_u8;
    }

    static void model(std::vector<uint8_t>& cntsu8, uint32_t const* in, uint32_t /*sum_of_values*/, size_t n)
    {
        auto counts_ptr = reinterpret_cast<ans_packed_counts*>(cntsu8.data());
        auto& counts = *counts_ptr;
        auto model_id = model_type::pick_model(in, n);
        uint32_t max_val = 0;
        for (size_t i = 0; i < n; i++) {
            uint8_t mag = ans_packed::magnitude(in[i] + 1);
            counts[model_id].counts[mag]++;
            max_val = std::max(max_val, in[i] + 1);
        }
        counts[model_id].max_value = std::max(max_val, counts[model_id].max_value);
    }

    static uint64_t encode_num_compact(const ans_packed::enc_model_compact* model, uint64_t state,
        uint32_t num, uint8_t*& out)
    {
        // (0) lookup quantities
        uint8_t mag = ans_packed::magnitude(num);
        uint64_t min_val = ans_packed::min_val_in_mag(mag);
        uint64_t vals_before = num - min_val;
        uint64_t freq = model->nfreq[mag];
        uint64_t base = model->base[mag] + (freq * vals_before);
        uint64_t SUB = model->SUB[mag];

        // (1) normalize
        while (state >= SUB) {
            ans_packed::output_unit<ans_packed::constants::OUTPUT_BASE_LOG2>(out, state);
        }

        // (2) transform state
        uint64_t next = ((state / freq) * model->M) + (state % freq) + base;
        return next;
    }

    static uint64_t encode_num(const ans_packed::enc_model* model, uint64_t state,
        uint32_t num, uint8_t*& out)
    {
        const auto& entry = model->table[num];
        uint32_t f = entry.freq;
        uint64_t b = entry.base;

        // (1) normalize
        while (state >= entry.SUB) {
            ans_packed::output_unit<ans_packed::constants::OUTPUT_BASE_LOG2>(out, state);
        }

        // (2) transform state
        uint64_t next = ((state / f) * model->M) + (state % f) + b;
        return next;
    }

    static void flush_state(uint64_t state, uint8_t*& out, size_t num_bytes)
    {
        for (size_t i = 0; i < num_bytes; i++) {
            uint8_t out_byte = state & 0xFF;
            out--;
            *out = out_byte;
            state >>= 8;
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

        // (2) reverse encode the block using the selected ANS model
        auto model_ptrs = reinterpret_cast<const uint64_t*>(enc_model_u8.data());
        size_t model_offset = model_ptrs[bh.model_id];
        uint64_t state = 0;
        auto out_ptr = tmp_out_buf.data() + tmp_out_buf.size() - 1;
        auto out_start = out_ptr;
        if (is_compact_model(model_offset)) {
            auto cur_model = reinterpret_cast<const ans_packed::enc_model_compact*>(enc_model_u8.data() + model_offset);
            for (size_t k = 0; k < n; k++) {
                uint32_t num = in[n - k - 1] + 1;
                state = encode_num_compact(cur_model, state, num, out_ptr);
            }
        } else {
            auto cur_model = reinterpret_cast<const ans_packed::enc_model*>(enc_model_u8.data() + model_offset);
            for (size_t k = 0; k < n; k++) {
                uint32_t num = in[n - k - 1] + 1;
                state = encode_num(cur_model, state, num, out_ptr);
            }
        }

        size_t enc_size = out_start - out_ptr;
        size_t u32s_written = enc_size / sizeof(uint32_t);

        // (3) write block header
        bh.final_state_bytes = ans_packed::state_bytes(state);
        bh.num_ans_u32s = u32s_written;
        model_type::write_block_header(bh, out);

        // (4) write the final state
        flush_state(state, out_ptr, bh.final_state_bytes);

        // (5) copy the ans output to the buffer
        size_t final_enc_size = out_start - out_ptr;
        out.insert(out.end(), out_ptr, out_ptr + final_enc_size);
    }

    static uint8_t find_mag(const ans_packed::dec_model_compact* model, uint64_t state_mod_M)
    {

        size_t i = 0;
        for (; i <= ans_packed::constants::MAX_MAG && model->base[i] == 0; i++) {
        }
        size_t last_change = i - 1;
        for (; i <= ans_packed::constants::MAX_MAG; i++) {
            if (model->base[i] > state_mod_M || model->base[i] == 0)
                return last_change;
            if (model->base[i] != model->base[last_change])
                last_change = i;
        }
        return ans_packed::constants::MAX_MAG;
    }

    static uint32_t decode_num_compact(const ans_packed::dec_model_compact* model, uint64_t& state, const uint8_t*& in, size_t& enc_size)
    {
        uint64_t state_mod_M = state & model->mask_M;
        uint8_t state_mag = find_mag(model, state_mod_M);
        uint32_t freq = model->nfreq[state_mag];
        uint64_t mag_offset = (state_mod_M - model->base[state_mag]);
        uint64_t offset = mag_offset % freq;
        uint64_t num_offset = mag_offset / freq;
        uint32_t num = ans_packed::min_val_in_mag(state_mag) + num_offset;
        state = freq * (state >> model->log2_M) + offset;
        while (enc_size && state < model->norm_lower_bound) {
            ans_packed::input_unit<ans_packed::constants::OUTPUT_BASE_LOG2>(in, state, enc_size);
        }

        return num;
    }

    static uint32_t decode_num(const ans_packed::dec_model* model, uint64_t& state, const uint8_t*& in, size_t& enc_size)
    {
        uint64_t state_mod_M = state & model->mask_M;
        const auto& entry = model->table[state_mod_M];
        uint32_t num = entry.sym;
        uint32_t f = entry.freq;
        state = f * (state >> model->log2_M) + entry.offset;
        while (enc_size && state < model->norm_lower_bound) {
            ans_packed::input_unit<ans_packed::constants::OUTPUT_BASE_LOG2>(in, state, enc_size);
        }

        return num;
    }

    static uint64_t init_decoder_state(const uint8_t*& in, uint8_t num_bytes)
    {
        uint64_t state = 0;
        for (size_t i = 0; i < num_bytes; i++) {
            uint8_t new_byte = *in++;
            state <<= 8;
            state = state + new_byte;
        }
        return state;
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

        uint64_t state = init_decoder_state(in, bh.final_state_bytes);

        if (is_compact_model(model_offset)) {
            auto cur_model = reinterpret_cast<const ans_packed::dec_model_compact*>(dec_model_u8 + model_offset);
            for (size_t k = 0; k < n; k++) {
                uint32_t dec_num = decode_num_compact(cur_model, state, in, enc_size);
                *out++ = dec_num - 1; // substract one as OT has 0s and our compactest num is 1
            }
        } else {
            auto cur_model = reinterpret_cast<const ans_packed::dec_model*>(dec_model_u8 + model_offset);
            for (size_t k = 0; k < n; k++) {
                uint32_t dec_num = decode_num(cur_model, state, in, enc_size);
                *out++ = dec_num - 1; // substract one as OT has 0s and our compactest num is 1
            }
        }
        return in;
    }
};
}
