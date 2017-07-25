#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "ans_decoding_stats.hpp"
#include "ans_msb_util.hpp"
#include "ans_util.hpp"
#include "block_codecs.hpp"

#define COMPACT_DEC_TABLE 1

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

    static std::vector<uint32_t> condense_models(std::vector<uint8_t>&)
    {
        std::vector<uint32_t> remapping(NUM_MODELS);
        for (size_t i = 0; i < NUM_MODELS; i++)
            remapping[i] = i;
        return remapping;
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

    static std::vector<uint32_t> condense_models(std::vector<uint8_t>&)
    {
        std::vector<uint32_t> remapping(NUM_MODELS);
        for (size_t i = 0; i < NUM_MODELS; i++)
            remapping[i] = i;
        return remapping;
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
        uint32_t model_id = (MAG2SEL[mag_90p] << 4) + MAG2SEL[mag_med];
        if (model_id == 0 && buf[n - 1] != 0)
            model_id++;
        return model_id;
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

    static std::vector<uint32_t> condense_models(std::vector<uint8_t>&)
    {
        std::vector<uint32_t> remapping(NUM_MODELS);
        for (size_t i = 0; i < NUM_MODELS; i++)
            remapping[i] = i;
        return remapping;
    }
};

struct msb_model_med90p_2d_merged {
    static const uint32_t NUM_MODELS = 16 * 16;
    static const uint32_t MAX_NUM_MODELS = 63;
    using ans_msb_counts_table = ans_msb::counts[NUM_MODELS];
    static uint32_t pick_model(uint32_t const* in, size_t n)
    {
        static const std::vector<uint32_t> MAG2SEL{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
            10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
        static std::vector<uint32_t> buf(ans::constants::BLOCK_SIZE);
        std::copy(in, in + n, buf.begin());
        std::sort(buf.begin(), buf.begin() + n);
        uint32_t mag_med = ans::magnitude(buf[n / 2] + 1);
        uint32_t mag_90p = ans::magnitude(buf[n * 0.9] + 1);
        uint32_t model_id = (MAG2SEL[mag_90p] << 4) + MAG2SEL[mag_med];
        if (model_id == 0 && buf[n - 1] != 0)
            model_id++;
        return model_id;
    }

    static void write_block_header(const msb_block_header& bh, std::vector<uint8_t>& out)
    {
        if (bh.model_id == 0) {
            out.push_back(0);
        } else {
            uint32_t header = (bh.model_id << 10) + ((bh.final_state_bytes - 1) << 7) + bh.num_ans_u32s;
            out.push_back(header >> 8);
            out.push_back(header & 0xFF);
        }
    }

    static void read_block_header(msb_block_header& bh, uint8_t const*& in)
    {
        uint8_t first_byte = *in++;
        if (first_byte == 0) {
            bh.model_id = 0;
            return;
        }
        uint16_t header = (uint16_t(first_byte) << 8) + *in++;
        bh.model_id = header >> 10;
        bh.final_state_bytes = ((header >> 7) & 0x7) + 1;
        bh.num_ans_u32s = header & 0x3F;
    }

    struct combine_cost {
        double H_a;
        double H_b;
    };

    static std::vector<uint32_t>
    condense_models(std::vector<uint8_t>& cntsu8)
    {
        auto counts_ptr = reinterpret_cast<ans_msb_counts_table*>(cntsu8.data());
        auto& counts = *counts_ptr;

        // (0) 0-out the 0th model. as it is special
        for (size_t k = 0; k <= ans_msb::constants::MAX_VAL; k++)
            counts[0][k] = 0;

        // (1) compute entropy for each individual model
        std::vector<uint32_t> remapping(NUM_MODELS, 0);
        std::vector<std::pair<double, uint64_t>> model_entropy(NUM_MODELS, { 0.0, 0 });
        size_t num_current_models = 0;
        for (size_t i = 0; i < NUM_MODELS; i++) {
            model_entropy[i] = ans_msb::compute_entropy(counts[i]);
            if (model_entropy[i].second != 0)
                num_current_models++;
        }

        // (2) merge the best models to reduce the total number of models
        std::vector<std::pair<uint64_t, uint64_t>> merge_ops;
        std::cout << "num_current_models = " << num_current_models << " target = " << MAX_NUM_MODELS << std::endl;
        while (num_current_models > MAX_NUM_MODELS) {
            // (a) compute all pairs
            auto min_pair = ans_msb::find_min_pair(counts, NUM_MODELS, model_entropy);

            // (b) remove one model and change the other
            ans_msb::merge_models(counts, model_entropy, min_pair.first, min_pair.second);

            // (c) keep track of what has moved where
            merge_ops.push_back(min_pair);

            num_current_models--;
        }

        // (3) create the remapping
        auto itr = merge_ops.rbegin();
        auto end = merge_ops.rend();
        while (itr != end) {
            auto op = *itr;
            auto from = op.first;
            auto to = op.second;
            std::cout << "process op from=" << from << " to=" << to << std::endl;
            if (remapping[to] != 0)
                to = remapping[to];
            remapping[from] = to;
            std::cout << "apply op from=" << from << " to=" << to << std::endl;
            ++itr;
        }

        // (4) to make the model ids small we have to remap
        //     to [1,63] and adjust the remapping again!
        std::vector<uint32_t> remapping_final(NUM_MODELS, 0);

        // (a) move the models into the right place
        size_t j = 1; // 0 is a reserved selector
        for (size_t i = 0; i < NUM_MODELS; i++) {
            if (remapping[i] == 0 && model_entropy[i].second != 0) { // actual model at this pos?
                remapping_final[i] = j;
                std::cout << "remap model " << i << " to " << j << std::endl;
                if (i != j) {
                    for (size_t k = 0; k <= ans_msb::constants::MAX_VAL; k++) {
                        counts[j][k] = counts[i][k];
                        counts[i][k] = 0;
                    }
                }
                j++;
            }
        }
        // (b) fix the redirects
        for (size_t i = 0; i < NUM_MODELS; i++) {
            if (remapping[i] != 0) { // model redirect
                remapping_final[i] = remapping_final[remapping[i]];
            }
        }

        return remapping_final;
    }
};

struct msb_model_medmax_2d_merged {
    static const uint32_t NUM_MODELS = 16 * 16;
    static const uint32_t MAX_NUM_MODELS = 63;
    using ans_msb_counts_table = ans_msb::counts[NUM_MODELS];
    static uint32_t pick_model(uint32_t const* in, size_t n)
    {
        static const std::vector<uint32_t> MAG2SEL{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
            10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
        static std::vector<uint32_t> buf(ans::constants::BLOCK_SIZE);
        memcpy(buf.data(), in, n * sizeof(uint32_t));
        //std::copy(in, in + n, buf.begin());
        std::sort(buf.begin(), buf.begin() + n);
        uint32_t mag_med = ans::magnitude(buf[n / 2] + 1);
        uint32_t mag_max = ans::magnitude(buf[n - 1] + 1);
        uint32_t model_id = (MAG2SEL[mag_max] << 4) + MAG2SEL[mag_med];
        if (model_id == 0 && buf[n - 1] != 0)
            model_id++;
        return model_id;
    }

    static void write_block_header(const msb_block_header& bh, std::vector<uint8_t>& out)
    {
        if (bh.model_id == 0) {
            out.push_back(0);
        } else {
            uint32_t header = (bh.model_id << 10) + ((bh.final_state_bytes - 1) << 7) + bh.num_ans_u32s;
            out.push_back(header >> 8);
            out.push_back(header & 0xFF);
        }
    }

    static void read_block_header(msb_block_header& bh, uint8_t const*& in)
    {
        uint8_t first_byte = *in++;
        if (first_byte == 0) {
            bh.model_id = 0;
            return;
        }
        uint16_t header = (uint16_t(first_byte) << 8) + *in++;
        bh.model_id = header >> 10;
        bh.final_state_bytes = ((header >> 7) & 0x7) + 1;
        bh.num_ans_u32s = header & 0x3F;
    }

    struct combine_cost {
        double H_a;
        double H_b;
    };

    static std::vector<uint32_t>
    condense_models(std::vector<uint8_t>& cntsu8)
    {
        auto counts_ptr = reinterpret_cast<ans_msb_counts_table*>(cntsu8.data());
        auto& counts = *counts_ptr;

        // (0) 0-out the 0th model. as it is special
        for (size_t k = 0; k <= ans_msb::constants::MAX_VAL; k++)
            counts[0][k] = 0;

        // (1) compute entropy for each individual model
        std::vector<uint32_t> remapping(NUM_MODELS, 0);
        std::vector<std::pair<double, uint64_t>> model_entropy(NUM_MODELS, { 0.0, 0 });
        size_t num_current_models = 0;
        for (size_t i = 0; i < NUM_MODELS; i++) {
            model_entropy[i] = ans_msb::compute_entropy(counts[i]);
            if (model_entropy[i].second != 0)
                num_current_models++;
        }

        // (2) merge the best models to reduce the total number of models
        std::vector<std::pair<uint64_t, uint64_t>> merge_ops;
        std::cout << "num_current_models = " << num_current_models << " target = " << MAX_NUM_MODELS << std::endl;
        while (num_current_models > MAX_NUM_MODELS) {
            // (a) compute all pairs
            auto min_pair = ans_msb::find_min_pair(counts, NUM_MODELS, model_entropy);

            // (b) remove one model and change the other
            ans_msb::merge_models(counts, model_entropy, min_pair.first, min_pair.second);

            // (c) keep track of what has moved where
            merge_ops.push_back(min_pair);

            num_current_models--;
        }

        // (3) create the remapping
        auto itr = merge_ops.rbegin();
        auto end = merge_ops.rend();
        while (itr != end) {
            auto op = *itr;
            auto from = op.first;
            auto to = op.second;
            std::cout << "process op from=" << from << " to=" << to << std::endl;
            if (remapping[to] != 0)
                to = remapping[to];
            remapping[from] = to;
            std::cout << "apply op from=" << from << " to=" << to << std::endl;
            ++itr;
        }

        // (4) to make the model ids small we have to remap
        //     to [1,63] and adjust the remapping again!
        std::vector<uint32_t> remapping_final(NUM_MODELS, 0);

        // (a) move the models into the right place
        size_t j = 1; // 0 is a reserved selector
        for (size_t i = 0; i < NUM_MODELS; i++) {
            if (remapping[i] == 0 && model_entropy[i].second != 0) { // actual model at this pos?
                remapping_final[i] = j;
                std::cout << "remap model " << i << " to " << j << std::endl;
                if (i != j) {
                    for (size_t k = 0; k <= ans_msb::constants::MAX_VAL; k++) {
                        counts[j][k] = counts[i][k];
                        counts[i][k] = 0;
                    }
                }
                j++;
            }
        }
        // (b) fix the redirects
        for (size_t i = 0; i < NUM_MODELS; i++) {
            if (remapping[i] != 0) { // model redirect
                remapping_final[i] = remapping_final[remapping[i]];
            }
        }

        return remapping_final;
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

    static std::vector<uint8_t> create_enc_model_from_counts(std::vector<uint8_t>& cntsu8)
    {
        // (0) condense models if necessary
        auto remap_models = model_type::condense_models(cntsu8);

        auto counts_ptr = reinterpret_cast<const ans_msb_counts_table*>(cntsu8.data());
        auto& counts = *counts_ptr;
        size_t pointers_to_models = NUM_MODELS * sizeof(uint32_t);
        size_t remap_table_size = NUM_MODELS * sizeof(uint32_t);
        std::vector<uint8_t> enc_models(pointers_to_models + remap_table_size, 0);

        // (1) copy the remap table
        auto remap_table_ptr = reinterpret_cast<uint32_t*>(enc_models.data()) + NUM_MODELS;
        std::cout << "REMAP [";
        for (size_t i = 0; i < NUM_MODELS; i++) {
            remap_table_ptr[i] = remap_models[i];
            std::cout << remap_models[i] << ",";
        }
        std::cout << "]" << std::endl;
        // (2) create the models
        for (size_t i = 0; i < NUM_MODELS; i++) {
            size_t model_offset = enc_models.size();
            // (1) create the model
            bool empty_model = ans_msb::create_enc_model(enc_models, counts[i]);
            auto model_offset_u32_ptr = reinterpret_cast<uint32_t*>(enc_models.data()) + i;
            // (1) store offset of model data in byte stream
            if (empty_model) {
                *model_offset_u32_ptr = 0;
            } else {
                *model_offset_u32_ptr = model_offset;
            }
        }

        return enc_models;
    }

    static std::vector<uint8_t> create_dec_models(const std::vector<uint8_t>& enc_models_u8)
    {
        auto enc_models = reinterpret_cast<const uint32_t*>(enc_models_u8.data());
        size_t pointers_to_models = NUM_MODELS * sizeof(uint32_t);
        std::vector<uint8_t> dec_models_u8(pointers_to_models, 0);
        for (size_t i = 0; i < NUM_MODELS; i++) {
            size_t dec_model_offset = dec_models_u8.size();
            if (enc_models[i] != 0) {
                size_t enc_model_offset = enc_models[i];
                auto enc_model_ptr = reinterpret_cast<const ans_msb::enc_model*>(enc_models_u8.data() + enc_model_offset);
                size_t model_size = dec_models_u8.size();
#ifdef COMPACT_DEC_TABLE
                ans_msb::create_dec_model_compact(dec_models_u8, *enc_model_ptr);
#else
                ans_msb::create_dec_model(dec_models_u8, *enc_model_ptr);
#endif
                model_size = dec_models_u8.size() - model_size;
                auto model_offset_u64_ptr = reinterpret_cast<uint32_t*>(dec_models_u8.data()) + i;
                *model_offset_u64_ptr = dec_model_offset;
            } else {
                auto model_offset_u64_ptr = reinterpret_cast<uint32_t*>(dec_models_u8.data()) + i;
                *model_offset_u64_ptr = 0;
            }
        }
        return dec_models_u8;
    }

    static void
    model(std::vector<uint8_t>& cntsu8, uint32_t const* in, uint32_t sum_of_values, size_t n)
    {
        auto counts_ptr = reinterpret_cast<ans_msb_counts_table*>(cntsu8.data());
        auto& counts = *counts_ptr;
        // exclude things we will not code using ANS
        if (sum_of_values != uint32_t(-1) && n <= ans_msb::constants::VBYTE_THRESHOLD) {
            return;
        }
        auto model_id = model_type::pick_model(in, n);
        for (size_t i = 0; i < n; i++) {
            auto msb_val = ans_msb::mapping_alistair(in[i] + 1);
            counts[model_id][msb_val]++;
        }
    }

    static void encode(uint32_t const* in, uint32_t sum_of_values,
        size_t n, std::vector<uint8_t>& out, const std::vector<uint8_t>& enc_model_u8)
    {
        if (sum_of_values == 0)
            return;
        static std::vector<uint8_t> tmp_out_buf(block_size * 8);
        static std::vector<uint8_t> exception_out_buf(block_size * 8);

        if (sum_of_values != uint32_t(-1) && n <= ans_msb::constants::VBYTE_THRESHOLD) {
            if (n == 1)
                return;
            size_t out_len = tmp_out_buf.size();
            TightVariableByte::encode(in, n, tmp_out_buf.data(), out_len);
            out.insert(out.end(), tmp_out_buf.data(), tmp_out_buf.data() + out_len);
            return;
        }

        // (1) determine and encode model id
        msb_block_header bh;
        bh.model_id = model_type::pick_model(in, n);

        // (2) remap the model id
        // std::cout << "model_id = " << bh.model_id << " remap to ";
        auto remap_table_ptr = reinterpret_cast<const uint32_t*>(enc_model_u8.data()) + NUM_MODELS;
        bh.model_id = remap_table_ptr[bh.model_id];
        // std::cout << bh.model_id << std::endl;

        if (bh.model_id == 0) { // all 1s. continue
            model_type::write_block_header(bh, out);
            return;
        }

        // (2) reverse encode the block using the selected ANS model
        auto model_ptrs = reinterpret_cast<const uint32_t*>(enc_model_u8.data());
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
        uint32_t sum_of_values, size_t n, uint8_t const* dec_model_u8)
    {
        if (sum_of_values == 0) {
            memset(out, 0, sizeof(uint32_t) * n);
            return in;
        }
        if (sum_of_values != uint32_t(-1) && n <= ans_msb::constants::VBYTE_THRESHOLD) {
            if (n == 1) {
                *out = sum_of_values;
                return in;
            }
            return TightVariableByte::decode(in, out, n);
        }

        msb_block_header bh;
        model_type::read_block_header(bh, in);

        // uniform block
        if (bh.model_id == 0) {
            memset(out, 0, sizeof(uint32_t) * n);
            return in;
        }
        size_t ans_enc_size = bh.num_ans_u32s * sizeof(uint32_t);
        auto model_ptrs = reinterpret_cast<const uint32_t*>(dec_model_u8);
        size_t model_offset = model_ptrs[bh.model_id];

        // (0) init the decoder
        uint64_t state = ans::init_decoder(in, bh.final_state_bytes);

        // (1) decode the ans parts and the exceptions at the same time
        auto except_ptr = in + ans_enc_size;
#ifdef COMPACT_DEC_TABLE
        auto cur_model = reinterpret_cast<const ans_msb::dec_model_small*>(dec_model_u8 + model_offset);
        auto cur_sym_table = reinterpret_cast<const ans_msb::dec_table_entry_small*>(cur_model->table_data);
        auto dec_sym_table_ptr = reinterpret_cast<const ans_msb::dec_table_entry_small_sym*>(cur_model->table_data + cur_model->M * sizeof(ans_msb::dec_table_entry_small));
#else
        auto cur_model = reinterpret_cast<const ans_msb::dec_model*>(dec_model_u8 + model_offset);
#endif
        size_t k = 0;
        for (; k < n; k += 4) {
#ifdef COMPACT_DEC_TABLE
            const auto& dec_entry1 = decode_num_compact(*cur_model, cur_sym_table, dec_sym_table_ptr, state, in, ans_enc_size);
#else
            const auto& dec_entry1 = decode_num(*cur_model, state, in, ans_enc_size);
#endif
#ifdef COMPACT_DEC_TABLE
            const auto& dec_entry2 = decode_num_compact(*cur_model, cur_sym_table, dec_sym_table_ptr, state, in, ans_enc_size);
#else
            const auto& dec_entry2 = decode_num(*cur_model, state, in, ans_enc_size);
#endif
#ifdef COMPACT_DEC_TABLE
            const auto& dec_entry3 = decode_num_compact(*cur_model, cur_sym_table, dec_sym_table_ptr, state, in, ans_enc_size);
#else
            const auto& dec_entry3 = decode_num(*cur_model, state, in, ans_enc_size);
#endif
#ifdef COMPACT_DEC_TABLE
            const auto& dec_entry4 = decode_num_compact(*cur_model, cur_sym_table, dec_sym_table_ptr, state, in, ans_enc_size);
#else
            const auto& dec_entry4 = decode_num(*cur_model, state, in, ans_enc_size);
#endif
            *out++ = ans_msb::undo_mapping(dec_entry1, except_ptr) - 1;
            *out++ = ans_msb::undo_mapping(dec_entry2, except_ptr) - 1;
            *out++ = ans_msb::undo_mapping(dec_entry3, except_ptr) - 1;
            *out++ = ans_msb::undo_mapping(dec_entry4, except_ptr) - 1;
        }

        for (; k < n; k++) {
#ifdef COMPACT_DEC_TABLE
            const auto& dec_entry = decode_num_compact(*cur_model, cur_sym_table, dec_sym_table_ptr, state, in, ans_enc_size);
#else
            const auto& dec_entry = decode_num(*cur_model, state, in, ans_enc_size);
#endif
            *out++ = ans_msb::undo_mapping(dec_entry, except_ptr) - 1;
        }
        return except_ptr;
    }
};
}
