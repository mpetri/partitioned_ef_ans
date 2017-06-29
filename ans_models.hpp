#pragma once

#include "block_codecs.hpp"

namespace constants {
const uint64_t ANS_START_STATE = 0;
const uint32_t OUTPUT_BASE = 256;
const uint8_t OUTPUT_BASE_LOG2 = 8;
const uint32_t MAX_SIGMA = 256;
const uint8_t MAX_MAG = 25;
const uint8_t NUM_MAGS = 16;
const std::array<uint8_t, NUM_MAGS> SEL2MAG{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12,
    14, 16, 19, 22, 25 };
const std::array<uint8_t, MAX_MAG + 1> MAG2SEL{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9,
    10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15 };
const uint64_t TOPFREQ = 1048576;
const uint64_t MAXSTACKSIZE = 10000000;
}

inline uint8_t ans_vbyte_size(uint64_t x)
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
    } else {
        return 9;
    }
    return 9;
}

inline uint64_t ans_vbyte_decode_u32(const uint8_t*& input, uint32_t& enc_size)
{
    uint64_t x = 0;
    uint64_t shift = 0;
    while (true) {
        uint8_t c = *input++;
        enc_size--;
        x += (uint32_t(c & 127) << shift);
        if (!(c & 128)) {
            return x;
        }
        shift += 7;
    }
    return x;
}

template <uint32_t i>
inline uint8_t ans_extract7bits(const uint32_t val)
{
    uint8_t v = static_cast<uint8_t>((val >> (7 * i)) & ((1ULL << 7) - 1));
    return v;
}

template <uint32_t i>
inline uint8_t ans_extract7bitsmaskless(const uint32_t val)
{
    uint8_t v = static_cast<uint8_t>((val >> (7 * i)));
    return v;
}

inline void ans_vbyte_freq_count(uint32_t x, uint64_t*& f)
{
    if (x < (1U << 7)) {
        f[x & 127]++;
    } else if (x < (1U << 14)) {
        f[ans_extract7bits<0>(x) | 128]++;
        f[ans_extract7bitsmaskless<1>(x) & 127]++;
    } else if (x < (1U << 21)) {
        f[ans_extract7bits<0>(x) | 128]++;
        f[ans_extract7bits<1>(x) | 128]++;
        f[ans_extract7bitsmaskless<2>(x) & 127]++;
    } else if (x < (1U << 28)) {
        f[ans_extract7bits<0>(x) | 128]++;
        f[ans_extract7bits<1>(x) | 128]++;
        f[ans_extract7bits<2>(x) | 128]++;
        f[ans_extract7bitsmaskless<3>(x) & 127]++;
    } else {
        f[ans_extract7bits<0>(x) | 128]++;
        f[ans_extract7bits<1>(x) | 128]++;
        f[ans_extract7bits<2>(x) | 128]++;
        f[ans_extract7bits<3>(x) | 128]++;
        f[ans_extract7bitsmaskless<4>(x) & 127]++;
    }
}

inline void ans_vbyte_encode_u32(uint8_t*& out, uint32_t x)
{
    if (x < (1ULL << 7)) {
        *out++ = static_cast<uint8_t>(x & 127);
    } else if (x < (1ULL << 14)) {
        *out++ = ans_extract7bits<0>(x) | 128;
        *out++ = ans_extract7bitsmaskless<1>(x) & 127;
    } else if (x < (1ULL << 21)) {
        *out++ = ans_extract7bits<0>(x) | 128;
        *out++ = ans_extract7bits<1>(x) | 128;
        *out++ = ans_extract7bitsmaskless<2>(x) & 127;
    } else if (x < (1ULL << 28)) {
        *out++ = ans_extract7bits<0>(x) | 128;
        *out++ = ans_extract7bits<1>(x) | 128;
        *out++ = ans_extract7bits<2>(x) | 128;
        *out++ = ans_extract7bitsmaskless<3>(x) & 127;
    } else {
        *out++ = ans_extract7bits<0>(x) | 128;
        *out++ = ans_extract7bits<1>(x) | 128;
        *out++ = ans_extract7bits<2>(x) | 128;
        *out++ = ans_extract7bits<3>(x) | 128;
        *out++ = ans_extract7bitsmaskless<4>(x) & 127;
    }
}

inline std::vector<uint8_t> ans_vbyte_encode(const uint32_t* in, uint32_t n)
{
    std::vector<uint8_t> out(n * 4);
    size_t j = 0;
    for (uint32_t i = 0; i < n; i++) {
        uint32_t x = in[i];
        if (x < (1ULL << 7)) {
            out[j++] = static_cast<uint8_t>(x & 127);
        } else if (x < (1ULL << 14)) {
            out[j++] = ans_extract7bits<0>(x) | 128;
            out[j++] = ans_extract7bitsmaskless<1>(x) & 127;
        } else if (x < (1ULL << 21)) {
            out[j++] = ans_extract7bits<0>(x) | 128;
            out[j++] = ans_extract7bits<1>(x) | 128;
            out[j++] = ans_extract7bitsmaskless<2>(x) & 127;
        } else if (x < (1ULL << 28)) {
            out[j++] = ans_extract7bits<0>(x) | 128;
            out[j++] = ans_extract7bits<1>(x) | 128;
            out[j++] = ans_extract7bits<2>(x) | 128;
            out[j++] = ans_extract7bitsmaskless<3>(x) & 127;
        } else {
            out[j++] = ans_extract7bits<0>(x) | 128;
            out[j++] = ans_extract7bits<1>(x) | 128;
            out[j++] = ans_extract7bits<2>(x) | 128;
            out[j++] = ans_extract7bits<3>(x) | 128;
            out[j++] = ans_extract7bitsmaskless<4>(x) & 127;
        }
    }
    out.resize(j);
    return std::move(out);
}

bool is_power_of_two(uint64_t x) { return ((x != 0) && !(x & (x - 1))); }

void ans_normalize_counts_power_of_two(const uint64_t* counts, size_t num, uint32_t* norm_counts, size_t target_power)
{
    for (size_t i = 0; i < num; i++)
        norm_counts[i] = counts[i];

    uint32_t n = 0;
    uint64_t initial_sum = 0;
    for (size_t i = 0; i < num; i++) {
        if (norm_counts[i] != 0) {
            n = i + 1;
            initial_sum += norm_counts[i];
        }
    }
    /* first phase in scaling process, distribute out the
       last bucket, assume it is the smallest n(s) area, scale
       the rest by the same amount */
    double fudge_factor = 0.95;
    uint64_t M = 0;
    while (true) {
        double C = double(target_power) / double(initial_sum);
        for (size_t i = 0; i < n; i++) {
            norm_counts[i] = fudge_factor * norm_counts[i] * C;
            if (counts[i] != 0 && norm_counts[i] < 1) {
                norm_counts[i] = 1;
            }
        }

        /* now, what does it all add up to? */
        M = 0;
        for (size_t m = 0; m < n; m++) {
            M += norm_counts[m];
        }
        /* do we have to try again? */
        if (M > target_power) {
            std::cout << "ans-normalization: first phase failed. M=" << M << " target_power=" << target_power << ". ";
            for (size_t i = 0; i < n; i++) {
                norm_counts[i] = counts[i];
            }
            fudge_factor -= 0.025;
            std::cout << "new fudge_factor=" << fudge_factor << std::endl;
        } else {
            break;
        }
    }
    /* fourth phase, round up to a power of two and then redistribute */
    uint64_t excess = target_power - M;

    /* flow that excess count backwards to the beginning of
       the selectors array, spreading it out across the buckets...
    */
    for (int64_t m = int64_t(n - 1); m >= 0; m--) {
        double ratio = double(excess) / double(M);
        uint64_t adder = ratio * norm_counts[m];
        if (adder > excess) {
            adder = excess;
        }
        excess -= adder;
        M -= norm_counts[m];
        norm_counts[m] += adder;
    }
    if (excess != 0) {
        norm_counts[0] += excess;
    }

    M = 0;
    for (size_t i = 0; i < n; i++) {
        M += norm_counts[i];
    }
    if (!is_power_of_two(M)) {
        std::cerr << "ERROR! not power of 2 after normalization = ("
                  << M << "," << target_power << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
}

namespace quasi_succinct {

struct ans_byte_enc_model {
    uint32_t M;
    uint32_t mask_M;
    uint32_t log2_M;
    uint32_t base[constants::MAX_SIGMA];
    uint32_t normalized_freqs[constants::MAX_SIGMA];
    uint32_t sym_upper_bound[constants::MAX_SIGMA];
};

struct dec_table_entry {
    uint16_t freq;
    uint16_t offset;
    uint32_t sym;
};

template <uint32_t t_frame_size>
struct ans_byte_dec_model {
    uint32_t mask_M;
    uint32_t log2_M;
    uint64_t norm_lower_bound;
    dec_table_entry table[t_frame_size];
};

template <uint32_t t_frame_size = 4096>
struct ans_vbyte_model {
    static const uint64_t block_size = 128;
    const uint32_t frame_size = t_frame_size;

    static std::vector<uint8_t> create_empty_counts()
    {
        size_t count_size = constants::MAX_SIGMA * sizeof(uint64_t);
        return std::vector<uint8_t>(count_size, 0);
    }

    static std::vector<uint8_t> create_enc_model_from_counts(const std::vector<uint8_t>& cntsu8)
    {
        auto counts = reinterpret_cast<const uint64_t*>(cntsu8.data());

        std::vector<uint8_t> enc_model(sizeof(ans_byte_enc_model));
        auto model = reinterpret_cast<ans_byte_enc_model*>(enc_model.data());
        model->M = t_frame_size;
        model->mask_M = t_frame_size - 1;
        model->log2_M = log2(t_frame_size);
        ans_normalize_counts_power_of_two(counts, constants::MAX_SIGMA, model->normalized_freqs, model->M);
        uint32_t cumsum = 0;
        for (size_t i = 0; i < constants::MAX_SIGMA; i++) {
            model->base[i] = cumsum;
            cumsum += model->normalized_freqs[i];
        }

        uint32_t norm_lower_bound = constants::OUTPUT_BASE * model->M;
        for (size_t j = 0; j < constants::MAX_SIGMA; j++) {
            model->sym_upper_bound[j]
                = ((norm_lower_bound / model->M) * constants::OUTPUT_BASE)
                * model->normalized_freqs[j];
        }

        return enc_model;
    }

    static std::vector<uint8_t> create_dec_model(const std::vector<uint8_t>& enc_model_u8)
    {
        auto enc_model = reinterpret_cast<const ans_byte_enc_model*>(enc_model_u8.data());
        std::vector<uint8_t> dec_model_u8(sizeof(ans_byte_dec_model<t_frame_size>));
        auto dec_model = reinterpret_cast<ans_byte_dec_model<t_frame_size>*>(dec_model_u8.data());
        dec_model->mask_M = enc_model->mask_M;
        dec_model->log2_M = enc_model->log2_M;
        dec_model->norm_lower_bound = constants::OUTPUT_BASE * enc_model->M;

        uint32_t base = 0;
        for (size_t j = 0; j < constants::MAX_SIGMA; j++) {
            uint16_t cur_freq = enc_model->normalized_freqs[j];
            for (size_t k = 0; k < cur_freq; k++) {
                dec_model->table[base + k].sym = j;
                dec_model->table[base + k].freq = cur_freq;
                dec_model->table[base + k].offset = k;
            }
            base += cur_freq;
        }

        return dec_model_u8;
    }

    static void model(std::vector<uint8_t>& cntsu8, uint32_t const* in, uint32_t /*sum_of_values*/, size_t n)
    {
        bool all_zero = true;
        for (size_t i = 0; i < n; i++) {
            if (in[i] != 0) {
                all_zero = false;
                break;
            }
        }
        if (all_zero)
            return;
        auto counts = reinterpret_cast<uint64_t*>(cntsu8.data());
        for (size_t i = 0; i < n; i++) {
            ans_vbyte_freq_count(in[i], counts);
        }
    }

    static uint32_t encode_sym(const ans_byte_enc_model* model, uint32_t state, uint8_t sym, uint8_t*& out8)
    {
        uint32_t f = model->normalized_freqs[sym];
        uint32_t b = model->base[sym];
        // (1) normalize
        uint32_t SUB = model->sym_upper_bound[sym];
        while (state >= SUB) {
            --out8;
            *out8 = (uint8_t)(state & 0xFF);
            state = state >> constants::OUTPUT_BASE_LOG2;
        }

        // (2) transform state
        uint32_t next = ((state / f) << model->log2_M) + (state % f) + b;
        assert(next == 0 || next > constants::OUTPUT_BASE * model->M);
        return next;
    }

    static void flush_state(uint32_t final_state, uint8_t*& out8)
    {
        auto vb_bytes = ans_vbyte_size(final_state);
        out8 -= vb_bytes;
        auto tmp = out8;
        ans_vbyte_encode_u32(tmp, final_state);
    }

    static void encode(uint32_t const* in, uint32_t /* sum_of_values */,
        size_t n, std::vector<uint8_t>& out, const std::vector<uint8_t>& enc_model_u8)
    {
        auto enc_model = reinterpret_cast<const ans_byte_enc_model*>(enc_model_u8.data());

        if (n != ans_vbyte_model::block_size) {
            static TightVariableByte vbyte_codec;
            std::vector<uint8_t> buf(2 * 4 * block_size);
            size_t out_len = buf.size();
            vbyte_codec.encode(in, n, buf.data(), out_len);
            out.insert(out.end(), buf.data(), buf.data() + out_len);
            return;
        }

        // (1) encode vbyte
        std::vector<uint8_t> tmp_vbyte_buf = ans_vbyte_encode(in, n);

        // (2) write num vbytes we will encode
        size_t num_vbytes = tmp_vbyte_buf.size();
        TightVariableByte::encode_single(num_vbytes - n, out);

        // (3) ans encode to tmp buf
        std::vector<uint8_t> buf(2 * 4 * n);
        auto tmp_out_ptr = buf.data() + buf.size() - 1;
        auto tmp_out_start = tmp_out_ptr;
        auto state = constants::ANS_START_STATE;
        auto encin = tmp_vbyte_buf.data() + num_vbytes - 1;
        for (size_t i = 0; i < num_vbytes; i++) {
            uint8_t sym = *encin--;
            state = encode_sym(enc_model, state, sym, tmp_out_ptr);
        }
        flush_state(state, tmp_out_ptr);
        // as we encoded in reverse order, we have to write in into a tmp
        // buf and then output the written bytes
        size_t enc_size = (tmp_out_start - tmp_out_ptr);
        TightVariableByte::encode_single(enc_size, out);

        // (4) copy to real out buf
        out.insert(out.end(), tmp_out_ptr, tmp_out_ptr + enc_size);
    }

    static uint32_t init_decoder(uint8_t const*& in, uint32_t& enc_size)
    {
        return ans_vbyte_decode_u32(in, enc_size);
    }

    static uint8_t const* decode(uint8_t const* in, uint32_t* out,
        uint32_t /* sum_of_values */, size_t n, uint8_t const* dec_model_u8)
    {
        auto model = reinterpret_cast<ans_byte_dec_model<t_frame_size> const*>(dec_model_u8);

        if (n != ans_vbyte_model::block_size) {
            static TightVariableByte vbyte_codec;
            return vbyte_codec.decode(in, out, n);
        }

        // (1) determine vbyte syms
        uint32_t num_vb;
        in = TightVariableByte::decode(in, &num_vb, 1);
        num_vb += n;

        // (2) read num encoded syms
        uint32_t enc_size = 0;
        in = TightVariableByte::decode(in, &enc_size, 1);
        uint32_t state = init_decoder(in, enc_size);

        uint8_t shift = 0;
        uint32_t cur_num = 0;
        for (uint32_t i = 0; i < num_vb; i++) {

            uint32_t state_mod_M = state & model->mask_M;
            const auto& entry = model->table[state_mod_M];

            // update state and renormalize
            state = entry.freq * (state >> model->log2_M) + entry.offset;
            while (enc_size && state < model->norm_lower_bound) {
                uint8_t new_byte = *in++;
                state = (state << constants::OUTPUT_BASE_LOG2) | uint32_t(new_byte);
                enc_size--;
            }

            uint8_t cur_sym = entry.sym;
            cur_num += (uint32_t(cur_sym & 127) << shift);
            if (!(cur_sym & 128)) {
                *out++ = cur_num;
                cur_num = 0;
                shift = 0;
            } else {
                shift += 7;
            }
        }

        return in;
    }
};
}
