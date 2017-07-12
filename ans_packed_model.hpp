#pragma once

#include "ans_util.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

using ans_packed_counts = ans_packed::mag_table[ans_packed::constants::NUM_MAGS];

namespace quasi_succinct {

struct ans_packed_model {
    static const uint64_t block_size = 128;

    static std::vector<uint8_t> create_empty_counts()
    {
        size_t count_size = sizeof(ans_packed_counts);
        return std::vector<uint8_t>(count_size, 0);
    }

    static ans_packed::mag_table* normalize_counts(const ans_packed::mag_table* table)
    {
        // print_mag_table(table, "initial_freqs");
        ans_packed::mag_table* nfreqs = new ans_packed::mag_table;
        nfreqs->max_value = table->max_value;
        uint64_t initial_sum = 0;
        for (size_t i = 0; i <= ans_packed::constants::MAX_MAG; i++) {
            initial_sum += table->counts[i] * ans_packed::uniq_vals_in_mag(i, table->max_value);
            nfreqs->counts[i] = table->counts[i];
        }

        uint8_t max_mag = 0;
        for (size_t i = 0; i <= ans_packed::constants::MAX_MAG; i++) {
            if (nfreqs->counts[i] != 0)
                max_mag = i;
        }
        /* first phase in scaling process, distribute out the
           last bucket, assume it is the smallest n(s) area, scale
           the rest by the same amount */
        auto bucket_size = ans_packed::uniq_vals_in_mag(max_mag, nfreqs->max_value);
        double C = 0.5 * bucket_size / nfreqs->counts[max_mag];
        for (size_t m = 0; m <= max_mag; m++) {
            bucket_size = ans_packed::uniq_vals_in_mag(m, nfreqs->max_value);
            nfreqs->counts[m] = 0.5 + nfreqs->counts[m] * C / bucket_size;
            if (table->counts[m] != 0 && nfreqs->counts[m] < 1) {
                nfreqs->counts[m] = 1;
            }
        }
        // print_mag_table(nfreqs, "first_phase");
        /* second step in scaling process, to make the first freq
           less than or equal to TOPFREQ
        */
        if (nfreqs->counts[0] > ans_packed::constants::TOPFREQ) {
            C = 1.0 * ans_packed::constants::TOPFREQ / nfreqs->counts[0];
            nfreqs->counts[0] = ans_packed::constants::TOPFREQ;
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
            M += nfreqs->counts[m] * ans_packed::uniq_vals_in_mag(m, nfreqs->max_value);
        }
        /* fourth phase, round up to a power of two and then redistribute */
        uint64_t target_power = ans_packed::next_power_of_two(M);
        uint64_t excess = target_power - M;
        /* flow that excess count backwards to the beginning of
           the selectors array, spreading it out across the buckets...
        */
        for (int8_t m = int8_t(max_mag); m >= 0; m--) {
            double ratio = 1.0 * excess / M;
            uint64_t adder = ratio * nfreqs->counts[m];
            excess -= ans_packed::uniq_vals_in_mag(m, nfreqs->max_value) * adder;
            M -= ans_packed::uniq_vals_in_mag(m, nfreqs->max_value) * nfreqs->counts[m];
            nfreqs->counts[m] += adder;
        }
        if (excess != 0) {
            nfreqs->counts[0] += excess;
        }
        // print_mag_table(nfreqs, "final_phase");

        M = 0;
        for (size_t i = 0; i <= max_mag; i++) {
            M += int64_t(nfreqs->counts[i] * ans_packed::uniq_vals_in_mag(i, nfreqs->max_value));
        }

        if (!ans_packed::is_power_of_two(M)) {
            fprintf(stderr, "ERROR! not power of 2 after normalization = %lu\n", M);
            exit(EXIT_FAILURE);
        }

        return nfreqs;
    }

    static bool create_enc_model(std::vector<uint8_t>& enc_models, const ans_packed::mag_table& table)
    {
        // (0) if all is 0 do nothing
        if (std::all_of(table.counts, table.counts + ans_packed::constants::MAX_MAG + 1,
                [](uint64_t i) { return i == 0; })) {
            return true;
        }

        // (1) normalize the counts
        auto norm_counts = normalize_counts(&table);

        // (2) create the encoding model
        size_t model_size = sizeof(ans_packed::enc_model) + (table.max_value + 1) * sizeof(ans_packed::mag_enc_table_entry);
        std::vector<uint8_t> new_model(model_size);
        auto model_ptr = reinterpret_cast<ans_packed::enc_model*>(new_model.data());
        auto& model = *model_ptr;

        // (2a) fill the tables
        uint64_t cumsum = 0;
        for (size_t i = 0; i <= ans_packed::constants::MAX_MAG; i++) {
            if (norm_counts->counts[i] == 0)
                continue;
            auto min_val = ans_packed::min_val_in_mag(i);
            auto max_val = ans_packed::max_val_in_mag(i, norm_counts->max_value);
            for (size_t j = min_val; j <= max_val; j++) {
                model.table[j].freq = norm_counts->counts[i];
                model.table[j].base = cumsum;
                cumsum += model.table[j].freq;
            }
        }

        model.M = cumsum;
        model.norm_lower_bound = ans_packed::constants::NORM_LOWER_BOUND;
        for (size_t j = 1; j < (norm_counts->max_value + 1); j++) {
            model.table[j].SUB = ((model.norm_lower_bound / model.M) * ans_packed::constants::OUTPUT_BASE)
                * model.table[j].freq;
        }
        model.mask_M = model.M - 1;
        model.log2_M = log2(model.M);
        model.max_value = norm_counts->max_value;

        enc_models.insert(enc_models.end(), new_model.begin(), new_model.end());
        delete norm_counts;
        return false;
    }

    static std::vector<uint8_t> create_enc_model_from_counts(const std::vector<uint8_t>& cntsu8)
    {
        auto counts_ptr = reinterpret_cast<const ans_packed_counts*>(cntsu8.data());
        const ans_packed_counts& counts = *counts_ptr;
        size_t pointers_to_models = ans_packed::constants::NUM_MAGS * sizeof(uint64_t);
        std::vector<uint8_t> enc_models(pointers_to_models, 0);

        for (size_t i = 0; i < ans_packed::constants::NUM_MAGS; i++) {
            // (1) store offset of model data in byte stream
            size_t model_offset = enc_models.size();
            // (2) create the model
            bool empty_model = create_enc_model(enc_models, counts[i]);
            auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(enc_models.data()) + i;
            if (empty_model) {
                *model_offset_u64_ptr = 0;
            } else {
                *model_offset_u64_ptr = model_offset;
            }
        }
        return enc_models;
    }

    static void create_dec_model(std::vector<uint8_t>& dec_models, const ans_packed::enc_model& enc_model)
    {
        // (1) determine model size
        size_t model_size = sizeof(ans_packed::dec_model) + (enc_model.M) * sizeof(ans_packed::mag_dec_table_entry);
        std::vector<uint8_t> new_model(model_size);
        auto model_ptr = reinterpret_cast<ans_packed::dec_model*>(new_model.data());
        auto& model = *model_ptr;

        // (2) create csum table for decoding
        model.M = enc_model.M;
        model.mask_M = model.M - 1;
        model.log2_M = log2(model.M);
        model.norm_lower_bound = enc_model.norm_lower_bound;
        size_t base = 0;
        for (size_t j = 1; j <= enc_model.max_value; j++) {
            auto cur_freq = enc_model.table[j].freq;
            for (size_t k = 0; k < cur_freq; k++) {
                model.table[base + k].sym = j;
                model.table[base + k].freq = cur_freq;
                model.table[base + k].offset = k;
            }
            base += cur_freq;
        }
        dec_models.insert(dec_models.end(), new_model.begin(), new_model.end());
    }

    static std::vector<uint8_t> create_dec_model(const std::vector<uint8_t>& enc_models_u8)
    {
        auto enc_models = reinterpret_cast<const uint64_t*>(enc_models_u8.data());
        size_t pointers_to_models = ans_packed::constants::NUM_MAGS * sizeof(uint64_t);
        std::vector<uint8_t> dec_models_u8(pointers_to_models, 0);

        for (size_t i = 0; i < ans_packed::constants::NUM_MAGS; i++) {
            size_t dec_model_offset = dec_models_u8.size();
            if (enc_models[i] != 0) {
                size_t enc_model_offset = enc_models[i];
                auto enc_model_ptr = reinterpret_cast<const ans_packed::enc_model*>(enc_models_u8.data() + enc_model_offset);
                const ans_packed::enc_model& enc_model = *enc_model_ptr;
                create_dec_model(dec_models_u8, enc_model);
                auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(dec_models_u8.data()) + i;
                *model_offset_u64_ptr = dec_model_offset;
            } else {
                auto model_offset_u64_ptr = reinterpret_cast<uint64_t*>(dec_models_u8.data()) + i;
                *model_offset_u64_ptr = 0;
            }
        }
        return dec_models_u8;
    }

    static uint8_t pick_model(uint32_t const* in, size_t n)
    {
        uint32_t max_val = 0;
        for (size_t i = 0; i < n; i++) {
            max_val = std::max(max_val, in[i] + 1);
        }
        uint8_t max_mag = ans_packed::magnitude(max_val);
        return ans_packed::constants::MAG2SEL[max_mag];
    }

    static void model(std::vector<uint8_t>& cntsu8, uint32_t const* in, uint32_t /*sum_of_values*/, size_t n)
    {
        auto counts_ptr = reinterpret_cast<ans_packed_counts*>(cntsu8.data());
        auto& counts = *counts_ptr;
        auto model_id = pick_model(in, n);
        uint32_t max_val = 0;
        for (size_t i = 0; i < n; i++) {
            uint8_t mag = ans_packed::magnitude(in[i] + 1);
            counts[model_id].counts[mag]++;
            max_val = std::max(max_val, in[i] + 1);
        }
        counts[model_id].max_value = std::max(max_val, counts[model_id].max_value);
    }

    static uint64_t encode_num(const ans_packed::enc_model* model, uint64_t state,
        uint32_t num, uint8_t*& out, bool start_block)
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
        auto model_id = pick_model(in, n);

        if (model_id == 0) { // all 1s. continue
            uint8_t packed = ans_packed::pack_two_4bit_nums(model_id, 0);
            out.push_back(packed);
            return;
        }

        static std::array<uint8_t, block_size * 8> tmp_out_buf;

        // (2) reverse encode the block using the selected ANS model
        auto model_ptrs = reinterpret_cast<const uint64_t*>(enc_model_u8.data());
        size_t model_offset = model_ptrs[model_id];
        auto cur_model = reinterpret_cast<const ans_packed::enc_model*>(enc_model_u8.data() + model_offset);
        uint64_t state = 0;
        auto out_ptr = tmp_out_buf.data() + tmp_out_buf.size() - 1;
        auto out_start = out_ptr;
        for (size_t k = 0; k < n; k++) {
            uint32_t num = in[n - k - 1] + 1;
            state = encode_num(cur_model, state, num, out_ptr, k == 0);
        }
        size_t enc_size = out_start - out_ptr;
        size_t u32s_written = enc_size / sizeof(uint32_t);

        // write model id and num final state bytes in 8bit
        size_t fsb = ans_packed::state_bytes(state);
        uint8_t packed = ans_packed::pack_two_4bit_nums(model_id, fsb);
        out.push_back(packed);

        // write the number of u32s of the output using vbyte
        // likely 1 byte only
        out.push_back(uint8_t(u32s_written));

        // (3) copy to real ANS out buf
        flush_state(state, out_ptr, fsb);
        size_t final_enc_size = out_start - out_ptr;
        out.insert(out.end(), out_ptr, out_ptr + final_enc_size);
    }

    static uint32_t decode_num(const ans_packed::dec_model* model, uint64_t& state, const uint8_t*& in, size_t& enc_size)
    {
        size_t state_start = state;
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
        uint8_t packed = *in++;
        auto model_and_fsb = ans_packed::unpack_two_4bit_nums(packed);
        uint8_t model_id = model_and_fsb.first;
        uint8_t fsb = model_and_fsb.second;

        // uniform block
        if (model_id == 0) {
            for (size_t i = 0; i < n; i++)
                out[i] = 0;
            return in;
        }

        uint8_t num_u32 = *in++;
        size_t enc_size = num_u32 * sizeof(uint32_t);

        auto model_ptrs = reinterpret_cast<const uint64_t*>(dec_model_u8);
        size_t model_offset = model_ptrs[model_id];
        auto cur_model = reinterpret_cast<const ans_packed::dec_model*>(dec_model_u8 + model_offset);

        uint64_t state = init_decoder_state(in, fsb);

        for (size_t k = 0; k < n; k++) {
            uint32_t dec_num = decode_num(cur_model, state, in, enc_size);
            *out++ = dec_num - 1; // substract one as OT has 0s and our smallest num is 1
        }
        return in;
    }
};
}
