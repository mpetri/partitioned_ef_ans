#pragma once

#include "block_codecs.hpp"

namespace ans_packed_constants {
const uint8_t MAX_MAG = 32;
const uint8_t NUM_MAGS = 16;
const std::array<uint8_t, NUM_MAGS> SEL2MAG{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12,
    14, 16, 19, 22, 32 };
const std::array<uint8_t, MAX_MAG + 1> MAG2SEL{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
    10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 };
const uint64_t TOPFREQ = 1048576;
const uint8_t OUTPUT_BASE_LOG2 = 32;
const uint64_t OUTPUT_BASE = 1ULL << OUTPUT_BASE_LOG2;
const uint64_t NORM_LOWER_BOUND = 1ULL << 24;
}

//#define ANS_DEBUG 1
static int ANS_DEBUG_IN = 1000;
static int ANS_DEBUG_OUT = 1;

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
void ans_packed_output_unit(uint8_t*& out, uint64_t& state)
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
void ans_packed_output_unit<8>(uint8_t*& out, uint64_t& state)
{
    --out;
    *out = (uint8_t)(state & 0xFF);
    state = state >> 8;
}

template <>
void ans_packed_output_unit<16>(uint8_t*& out, uint64_t& state)
{
    out -= 2;
    uint16_t* out16 = reinterpret_cast<uint16_t*>(out);
    *out16 = (uint16_t)(state & 0xFFFF);
    state = state >> 16;
}

template <>
void ans_packed_output_unit<32>(uint8_t*& out, uint64_t& state)
{
    out -= 4;
    uint32_t* out32 = reinterpret_cast<uint32_t*>(out);
    *out32 = (uint32_t)(state & 0xFFFFFFFF);
    state = state >> 32;
}

template <uint8_t t_width>
void ans_packed_input_unit(const uint8_t*& in, uint64_t& state, size_t& enc_size)
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
void ans_packed_input_unit<8>(const uint8_t*& in, uint64_t& state, size_t& enc_size)
{
    uint8_t new_byte = *in++;
    state = (state << 8) | uint64_t(new_byte);
    enc_size--;
}

template <>
void ans_packed_input_unit<16>(const uint8_t*& in, uint64_t& state, size_t& enc_size)
{
    const uint16_t* in16 = reinterpret_cast<const uint16_t*>(in);
    uint64_t new_unit = *in16;
    state = (state << 16) | new_unit;
    in += 2;
    enc_size -= 2;
}

template <>
void ans_packed_input_unit<32>(const uint8_t*& in, uint64_t& state, size_t& enc_size)
{
    const uint32_t* in32 = reinterpret_cast<const uint32_t*>(in);
    uint64_t new_unit = *in32;
    state = (state << 32) | new_unit;
    in += 4;
    enc_size -= 4;
}

inline uint8_t ans_packed_vb_size(uint64_t x)
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

inline uint64_t ans_packed_vbyte_decode_u64(const uint8_t*& input)
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
inline uint8_t ans_packed_extract7bits(const uint64_t val)
{
    uint8_t v = static_cast<uint8_t>((val >> (7 * i)) & ((1ULL << 7) - 1));
    return v;
}

template <uint32_t i>
inline uint8_t ans_packed_extract7bitsmaskless(const uint64_t val)
{
    uint8_t v = static_cast<uint8_t>((val >> (7 * i)));
    return v;
}

inline void ans_packed_vbyte_encode_u64(uint8_t*& out, uint64_t x)
{
    if (x < (1ULL << 7)) {
        *out++ = static_cast<uint8_t>(x & 127);
    } else if (x < (1ULL << 14)) {
        *out++ = ans_packed_extract7bits<0>(x) | 128;
        *out++ = ans_packed_extract7bitsmaskless<1>(x) & 127;
    } else if (x < (1ULL << 21)) {
        *out++ = ans_packed_extract7bits<0>(x) | 128;
        *out++ = ans_packed_extract7bits<1>(x) | 128;
        *out++ = ans_packed_extract7bitsmaskless<2>(x) & 127;
    } else if (x < (1ULL << 28)) {
        *out++ = ans_packed_extract7bits<0>(x) | 128;
        *out++ = ans_packed_extract7bits<1>(x) | 128;
        *out++ = ans_packed_extract7bits<2>(x) | 128;
        *out++ = ans_packed_extract7bitsmaskless<3>(x) & 127;
    } else if (x < (1ULL << 35)) {
        *out++ = ans_packed_extract7bits<0>(x) | 128;
        *out++ = ans_packed_extract7bits<1>(x) | 128;
        *out++ = ans_packed_extract7bits<2>(x) | 128;
        *out++ = ans_packed_extract7bits<3>(x) | 128;
        *out++ = ans_packed_extract7bitsmaskless<4>(x) & 127;
    } else if (x < (1ULL << 42)) {
        *out++ = ans_packed_extract7bits<0>(x) | 128;
        *out++ = ans_packed_extract7bits<1>(x) | 128;
        *out++ = ans_packed_extract7bits<2>(x) | 128;
        *out++ = ans_packed_extract7bits<3>(x) | 128;
        *out++ = ans_packed_extract7bits<4>(x) | 128;
        *out++ = ans_packed_extract7bitsmaskless<5>(x) & 127;
    } else if (x < (1ULL << 49)) {
        *out++ = ans_packed_extract7bits<0>(x) | 128;
        *out++ = ans_packed_extract7bits<1>(x) | 128;
        *out++ = ans_packed_extract7bits<2>(x) | 128;
        *out++ = ans_packed_extract7bits<3>(x) | 128;
        *out++ = ans_packed_extract7bits<4>(x) | 128;
        *out++ = ans_packed_extract7bits<5>(x) | 128;
        *out++ = ans_packed_extract7bitsmaskless<6>(x) & 127;
    } else if (x < (1ULL << 56)) {
        *out++ = ans_packed_extract7bits<0>(x) | 128;
        *out++ = ans_packed_extract7bits<1>(x) | 128;
        *out++ = ans_packed_extract7bits<2>(x) | 128;
        *out++ = ans_packed_extract7bits<3>(x) | 128;
        *out++ = ans_packed_extract7bits<4>(x) | 128;
        *out++ = ans_packed_extract7bits<5>(x) | 128;
        *out++ = ans_packed_extract7bits<6>(x) | 128;
        *out++ = ans_packed_extract7bitsmaskless<7>(x) & 127;
    } else {
        *out++ = ans_packed_extract7bits<0>(x) | 128;
        *out++ = ans_packed_extract7bits<1>(x) | 128;
        *out++ = ans_packed_extract7bits<2>(x) | 128;
        *out++ = ans_packed_extract7bits<3>(x) | 128;
        *out++ = ans_packed_extract7bits<4>(x) | 128;
        *out++ = ans_packed_extract7bits<5>(x) | 128;
        *out++ = ans_packed_extract7bits<6>(x) | 128;
        *out++ = ans_packed_extract7bits<7>(x) | 128;
        *out++ = ans_packed_extract7bitsmaskless<8>(x) & 127;
    }
}

uint8_t ans_packed_magnitude(uint32_t x)
{
    uint64_t y = x;
    if (x == 1)
        return 0;
    uint32_t res = 63 - __builtin_clzll(y);
    if ((1ULL << res) == y)
        return res;
    return res + 1;
}

uint32_t ans_max_val_in_mag(uint8_t mag, uint32_t max_val = 0)
{
    uint32_t maxv = 1;
    if (mag != 0)
        maxv = (1ULL << (mag));
    if (maxv > max_val)
        maxv = max_val;
    return maxv;
}

uint32_t ans_min_val_in_mag(uint8_t mag)
{
    if (mag == 0)
        return 1;
    return (1ULL << (mag - 1)) + 1;
}

uint32_t ans_uniq_vals_in_mag(uint8_t mag, uint32_t max_val = 0)
{
    return ans_max_val_in_mag(mag, max_val) - ans_min_val_in_mag(mag) + 1;
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

struct mag_enc_table_entry {
    uint32_t freq;
    uint64_t base;
    uint64_t SUB;
};

struct mag_dec_table_entry {
    uint32_t freq;
    uint64_t offset;
    uint32_t sym;
};

struct ans_packed_enc_model {
    uint64_t M = 0; // frame size
    uint8_t log2_M = 0;
    uint64_t mask_M = 0;
    uint64_t norm_lower_bound = 0;
    uint32_t max_value = 0;
    mag_enc_table_entry table[0];
};

struct ans_packed_mag_table {
    uint32_t max_value;
    uint64_t counts[ans_packed_constants::MAX_MAG + 1];
};

void print_mag_table(const ans_packed_mag_table* tb, std::string name)
{
    std::cout << name << " max_value = " << tb->max_value << " COUNTS = (";
    for (size_t i = 0; i <= ans_packed_constants::MAX_MAG; i++) {
        std::cout << "<" << i << "," << tb->counts[i] << ">";
    }
    std::cout << ")" << std::endl;
}

struct ans_packed_dec_model {
    uint64_t M = 0; // frame size
    uint8_t log2_M = 0;
    uint64_t mask_M = 0;
    uint64_t norm_lower_bound = 0;
    mag_dec_table_entry table[0];
};

using ans_packed_counts = ans_packed_mag_table[ans_packed_constants::NUM_MAGS];

namespace quasi_succinct {

struct ans_packed_model {
    static const uint64_t block_size = 128;

    static std::vector<uint8_t> create_empty_counts()
    {
        size_t count_size = sizeof(ans_packed_counts);
        return std::vector<uint8_t>(count_size, 0);
    }

    static ans_packed_mag_table* normalize_counts(const ans_packed_mag_table* table)
    {
        // print_mag_table(table, "initial_freqs");
        ans_packed_mag_table* nfreqs = new ans_packed_mag_table;
        nfreqs->max_value = table->max_value;
        uint64_t initial_sum = 0;
        for (size_t i = 0; i <= ans_packed_constants::MAX_MAG; i++) {
            initial_sum += table->counts[i] * ans_uniq_vals_in_mag(i, table->max_value);
            nfreqs->counts[i] = table->counts[i];
        }

        uint8_t max_mag = 0;
        for (size_t i = 0; i <= ans_packed_constants::MAX_MAG; i++) {
            if (nfreqs->counts[i] != 0)
                max_mag = i;
        }
        /* first phase in scaling process, distribute out the
           last bucket, assume it is the smallest n(s) area, scale
           the rest by the same amount */
        auto bucket_size = ans_uniq_vals_in_mag(max_mag, nfreqs->max_value);
        double C = 0.5 * bucket_size / nfreqs->counts[max_mag];
        for (size_t m = 0; m <= max_mag; m++) {
            bucket_size = ans_uniq_vals_in_mag(m, nfreqs->max_value);
            nfreqs->counts[m] = 0.5 + nfreqs->counts[m] * C / bucket_size;
            if (table->counts[m] != 0 && nfreqs->counts[m] < 1) {
                nfreqs->counts[m] = 1;
            }
        }
        // print_mag_table(nfreqs, "first_phase");
        /* second step in scaling process, to make the first freq
           less than or equal to TOPFREQ
        */
        if (nfreqs->counts[0] > ans_packed_constants::TOPFREQ) {
            C = 1.0 * ans_packed_constants::TOPFREQ / nfreqs->counts[0];
            nfreqs->counts[0] = ans_packed_constants::TOPFREQ;
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
            M += nfreqs->counts[m] * ans_uniq_vals_in_mag(m, nfreqs->max_value);
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
            excess -= ans_uniq_vals_in_mag(m, nfreqs->max_value) * adder;
            M -= ans_uniq_vals_in_mag(m, nfreqs->max_value) * nfreqs->counts[m];
            nfreqs->counts[m] += adder;
        }
        if (excess != 0) {
            nfreqs->counts[0] += excess;
        }
        // print_mag_table(nfreqs, "final_phase");

        M = 0;
        for (size_t i = 0; i <= max_mag; i++) {
            M += int64_t(nfreqs->counts[i] * ans_uniq_vals_in_mag(i, nfreqs->max_value));
        }

        if (!is_power_of_two(M)) {
            fprintf(stderr, "ERROR! not power of 2 after normalization = %lu\n", M);
            exit(EXIT_FAILURE);
        }

        return nfreqs;
    }

    static bool create_enc_model(std::vector<uint8_t>& enc_models, const ans_packed_mag_table& table)
    {
        // (0) if all is 0 do nothing
        if (std::all_of(table.counts, table.counts + ans_packed_constants::MAX_MAG + 1,
                [](uint64_t i) { return i == 0; })) {
            return true;
        }

        // (1) normalize the counts
        auto norm_counts = normalize_counts(&table);

        // (2) create the encoding model
        size_t model_size = sizeof(ans_packed_enc_model) + (table.max_value + 1) * sizeof(mag_enc_table_entry);
        std::vector<uint8_t> new_model(model_size);
        auto model_ptr = reinterpret_cast<ans_packed_enc_model*>(new_model.data());
        auto& model = *model_ptr;

        // (2a) fill the tables
        uint64_t cumsum = 0;
        for (size_t i = 0; i <= ans_packed_constants::MAX_MAG; i++) {
            if (norm_counts->counts[i] == 0)
                continue;
            auto min_val = ans_min_val_in_mag(i);
            auto max_val = ans_max_val_in_mag(i, norm_counts->max_value);
            for (size_t j = min_val; j <= max_val; j++) {
                model.table[j].freq = norm_counts->counts[i];
                model.table[j].base = cumsum;
                cumsum += model.table[j].freq;
            }
        }

        model.M = cumsum;
        model.norm_lower_bound = ans_packed_constants::NORM_LOWER_BOUND;
        for (size_t j = 1; j < (norm_counts->max_value + 1); j++) {
            model.table[j].SUB = ((model.norm_lower_bound / model.M) * ans_packed_constants::OUTPUT_BASE)
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
        size_t pointers_to_models = ans_packed_constants::NUM_MAGS * sizeof(uint64_t);
        std::vector<uint8_t> enc_models(pointers_to_models, 0);

        for (size_t i = 0; i < ans_packed_constants::NUM_MAGS; i++) {
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

    static void create_dec_model(std::vector<uint8_t>& dec_models, const ans_packed_enc_model& enc_model)
    {
        // (1) determine model size
        size_t model_size = sizeof(ans_packed_dec_model) + (enc_model.M) * sizeof(mag_dec_table_entry);
        std::vector<uint8_t> new_model(model_size);
        auto model_ptr = reinterpret_cast<ans_packed_dec_model*>(new_model.data());
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
        size_t pointers_to_models = ans_packed_constants::NUM_MAGS * sizeof(uint64_t);
        std::vector<uint8_t> dec_models_u8(pointers_to_models, 0);

        for (size_t i = 0; i < ans_packed_constants::NUM_MAGS; i++) {
            size_t dec_model_offset = dec_models_u8.size();
            if (enc_models[i] != 0) {
                size_t enc_model_offset = enc_models[i];
                auto enc_model_ptr = reinterpret_cast<const ans_packed_enc_model*>(enc_models_u8.data() + enc_model_offset);
                const ans_packed_enc_model& enc_model = *enc_model_ptr;
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
        uint8_t max_mag = ans_packed_magnitude(max_val);
        return ans_packed_constants::MAG2SEL[max_mag];
    }

    static void model(std::vector<uint8_t>& cntsu8, uint32_t const* in, uint32_t /*sum_of_values*/, size_t n)
    {
        auto counts_ptr = reinterpret_cast<ans_packed_counts*>(cntsu8.data());
        auto& counts = *counts_ptr;
        auto model_id = pick_model(in, n);
        uint32_t max_val = 0;
        for (size_t i = 0; i < n; i++) {
            uint8_t mag = ans_packed_magnitude(in[i] + 1);
            counts[model_id].counts[mag]++;
            max_val = std::max(max_val, in[i] + 1);
        }
        counts[model_id].max_value = std::max(max_val, counts[model_id].max_value);
    }

    static uint64_t encode_num(const ans_packed_enc_model* model, uint64_t state,
        uint32_t num, uint8_t*& out, bool start_block)
    {
        const auto& entry = model->table[num];
        uint32_t f = entry.freq;
        uint64_t b = entry.base;

#ifdef ANS_DEBUG
        uint64_t start_state = state;
        static double total_expected_bits = 0;
        static double total_actual_bits = 0;
        static double total_emitted_bits = 0;
        if (start_block) {
            total_expected_bits = 0;
            total_actual_bits = 0;
            total_emitted_bits = 0;
        }
        double expected_increase = double(model->M) / double(f);
        double expected_bits = 0;
        if (expected_increase != 0)
            expected_bits = log2(expected_increase);

        total_expected_bits += expected_bits;
        uint32_t bits_emitted = 0;
#endif

        // (1) normalize
        while (state >= entry.SUB) {
            ans_packed_output_unit<ans_packed_constants::OUTPUT_BASE_LOG2>(out, state);
#ifdef ANS_DEBUG
            bits_emitted += ans_packed_constants::OUTPUT_BASE_LOG2;
#endif
        }

#ifdef ANS_DEBUG
        total_emitted_bits += bits_emitted;
        size_t state_after_norm = state;
#endif

        // (2) transform state
        uint64_t next = ((state / f) * model->M) + (state % f) + b;
#ifdef ANS_DEBUG
        double actual_bits = log2(next);
        if (state_after_norm != 0)
            actual_bits -= log2(state_after_norm);
        total_actual_bits += actual_bits;

        if (ANS_DEBUG_IN > 0) {
            std::cout << "START ANS_ENCODE(num=" << num << ",freq="
                      << entry.freq << ",base=" << b << ",SUB="
                      << entry.SUB << ",M=" << model->M << ",L=" << model->norm_lower_bound << ")\n";
            std::cout << "\tstart_state = " << start_state << "\n";
            std::cout << "\texpected_state_increase = " << expected_increase << "\n";
            std::cout << "\texpected_bits_increase = " << expected_bits << "\n";
            std::cout << "\tstate_after_normalization = " << state_after_norm << "\n";
            std::cout << "\tbits_emitted = " << bits_emitted << "\n";
            std::cout << "\tactual_state_increase = " << next - state_after_norm << "\n";
            std::cout << "\tnew_state = " << next << "\n";
            std::cout << "\tactual_bits_increase = " << actual_bits << "\n";
            std::cout << "\ttotal_expected_bits = " << total_expected_bits << "\n";
            std::cout << "\ttotal_actual_bits = " << total_actual_bits << "\n";
            std::cout << "\ttotal_emitted_bits = " << total_emitted_bits << "\n";
            std::cout << "STOP ANS_ENCODE(num=" << num << ")\n";
        }
#endif

        return next;
    }

    static void flush_state(uint64_t state, uint8_t*& out, size_t num_bytes)
    {
#ifdef ANS_DEBUG
        if (ANS_DEBUG_IN > 0) {
            std::cout << "START FLUSH_STATE(state=" << state << ")\n";
            std::cout << "\texpected_final_state_bits = " << log2(state) << "\n";
            std::cout << "\tactual_final_state_bits = " << num_bytes * 8 << "\n";
        }
#endif
        for (size_t i = 0; i < num_bytes; i++) {
            uint8_t out_byte = state & 0xFF;
            out--;
            *out = out_byte;
            state >>= 8;
        }
#ifdef ANS_DEBUG
        if (ANS_DEBUG_IN > 0) {
            std::cout << "STOP FLUSH_STATE(state=" << state << ")\n";
        }
#endif
    }

    static void encode(uint32_t const* in, uint32_t /*sum_of_values*/,
        size_t n, std::vector<uint8_t>& out, const std::vector<uint8_t>& enc_model_u8)
    {
        // (1) determine and encode model id
        auto model_id = pick_model(in, n);

        if (model_id == 0) { // all 1s. continue
            uint8_t packed = pack_two_4bit_nums(model_id, 0);
            out.push_back(packed);
            return;
        }

        static std::array<uint8_t, block_size * 8> tmp_out_buf;

        // (2) reverse encode the block using the selected ANS model
        auto model_ptrs = reinterpret_cast<const uint64_t*>(enc_model_u8.data());
        size_t model_offset = model_ptrs[model_id];
        auto cur_model = reinterpret_cast<const ans_packed_enc_model*>(enc_model_u8.data() + model_offset);
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
        size_t fsb = state_bytes(state);
        uint8_t packed = pack_two_4bit_nums(model_id, fsb);
        out.push_back(packed);

        // write the number of u32s of the output using vbyte
        // likely 1 byte only
        TightVariableByte::encode_single(u32s_written, out);

        // (3) copy to real ANS out buf
        flush_state(state, out_ptr, fsb);
        size_t final_enc_size = out_start - out_ptr;
        out.insert(out.end(), out_ptr, out_ptr + final_enc_size);
        if (ANS_DEBUG_IN)
            ANS_DEBUG_IN--;
    }

    static uint32_t decode_num(const ans_packed_dec_model* model, uint64_t& state, const uint8_t*& in, size_t& enc_size)
    {
        size_t state_start = state;
        uint64_t state_mod_M = state & model->mask_M;
        const auto& entry = model->table[state_mod_M];
        uint32_t num = entry.sym;
        uint32_t f = entry.freq;
        state = f * (state >> model->log2_M) + entry.offset;
        while (enc_size && state < model->norm_lower_bound) {
            ans_packed_input_unit<ans_packed_constants::OUTPUT_BASE_LOG2>(in, state, enc_size);
        }

#ifdef ANS_DEBUG
        if (ANS_DEBUG_OUT == 1) {
            std::cout << "START ANS_DECODE_NUM()" << std::endl;
            std::cout << "\tstate_start = " << state_start << std::endl;
            std::cout << "\tstate = " << state << std::endl;
            std::cout << "\tnorm_lower_bound = " << model->norm_lower_bound << std::endl;
            std::cout << "\tenc_size = " << enc_size << std::endl;
            std::cout << "STOP ANS_DECODE_NUM()" << std::endl;
        }
#endif
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
        auto model_and_fsb = unpack_two_4bit_nums(packed);
        uint8_t model_id = model_and_fsb.first;
        uint8_t fsb = model_and_fsb.second;

        // uniform block
        if (model_id == 0) {
            for (size_t i = 0; i < n; i++)
                out[i] = 0;
            return in;
        }

        uint32_t num_u32;
        in = TightVariableByte::decode(in, &num_u32, 1);
        size_t enc_size = num_u32 * sizeof(uint32_t);

        auto model_ptrs = reinterpret_cast<const uint64_t*>(dec_model_u8);
        size_t model_offset = model_ptrs[model_id];
        auto cur_model = reinterpret_cast<const ans_packed_dec_model*>(dec_model_u8 + model_offset);

        uint64_t state = init_decoder_state(in, fsb);

#ifdef ANS_DEBUG
        if (ANS_DEBUG_OUT == 1) {
            std::cout << "START ANS_DECODE\n";
            std::cout << "model_id = " << (int)model_id << std::endl;
            std::cout << "fsb = " << (int)fsb << std::endl;
            std::cout << "enc_size = " << enc_size << std::endl;
            std::cout << "initial_state = " << state << std::endl;
        }
#endif

        for (size_t k = 0; k < n; k++) {
            uint32_t dec_num = decode_num(cur_model, state, in, enc_size);
            *out++ = dec_num - 1; // substract one as OT has 0s and our smallest num is 1
        }

        ANS_DEBUG_OUT = 0;
        return in;
    }
};
}
