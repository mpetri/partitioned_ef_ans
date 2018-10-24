#pragma once

#include "block_codecs.hpp"
#include "succinct/util.hpp"
#include "util.hpp"
#include "ans_util.hpp"

namespace quasi_succinct {

struct ans_block_size_stats {
    ans_block_size_stats()
    {
        last_overhead_bytes = 0;
        small_overhead_bytes = 0;
        full_overhead_bytes = 0;
        total_overhead_bytes = 0;
        total_doc_bytes = 0;
        total_freq_bytes = 0;
        full_doc_bytes = 0;
        full_freq_bytes = 0;
        small_list_doc_bytes = 0;
        small_list_freq_bytes = 0;
        last_nonfull_doc_bytes = 0;
        last_nonfull_freq_bytes = 0;
        total_postings = 0;
        small_list_postings = 0;
        last_nonfull_postings = 0;
        full_block_postings = 0;
        for (size_t i = 0; i < ans::constants::BLOCK_SIZE; i++) {
            small_list_doc_bytesA[i] = 0;
            small_list_freq_bytesA[i] = 0;
            small_list_doc_postingsA[i] = 0;
            small_list_freq_postingsA[i] = 0;
        }
    }

    uint64_t last_overhead_bytes;
    uint64_t small_overhead_bytes;
    uint64_t full_overhead_bytes;
    uint64_t total_overhead_bytes;
    uint64_t total_doc_bytes;
    uint64_t total_freq_bytes;
    uint64_t full_doc_bytes;
    uint64_t full_freq_bytes;
    uint64_t small_list_doc_bytes;
    uint64_t small_list_freq_bytes;
    uint64_t last_nonfull_doc_bytes;
    uint64_t last_nonfull_freq_bytes;
    uint64_t total_postings;
    uint64_t small_list_postings;
    uint64_t last_nonfull_postings;
    uint64_t full_block_postings;
    uint64_t small_list_doc_bytesA[ans::constants::BLOCK_SIZE];
    uint64_t small_list_freq_bytesA[ans::constants::BLOCK_SIZE];
    uint64_t small_list_doc_postingsA[ans::constants::BLOCK_SIZE];
    uint64_t small_list_freq_postingsA[ans::constants::BLOCK_SIZE];

    ans_block_size_stats& operator+=(const ans_block_size_stats& rhs)
    {
        this->last_overhead_bytes += rhs.last_overhead_bytes;
        this->small_overhead_bytes += rhs.small_overhead_bytes;
        this->full_overhead_bytes += rhs.full_overhead_bytes;
        this->total_overhead_bytes += rhs.total_overhead_bytes;
        this->total_postings += rhs.total_postings;
        this->total_doc_bytes += rhs.total_doc_bytes;
        this->total_freq_bytes += rhs.total_freq_bytes;
        this->full_block_postings += rhs.full_block_postings;
        this->full_doc_bytes += rhs.full_doc_bytes;
        this->full_freq_bytes += rhs.full_freq_bytes;
        this->small_list_postings += rhs.small_list_postings;
        this->small_list_doc_bytes += rhs.small_list_doc_bytes;
        this->small_list_freq_bytes += rhs.small_list_freq_bytes;
        this->last_nonfull_doc_bytes += rhs.last_nonfull_doc_bytes;
        this->last_nonfull_freq_bytes += rhs.last_nonfull_freq_bytes;
        this->last_nonfull_postings += rhs.last_nonfull_postings;
        for (size_t i = 0; i < ans::constants::BLOCK_SIZE; i++) {
            this->small_list_doc_bytesA[i] += rhs.small_list_doc_bytesA[i];
            this->small_list_freq_bytesA[i] += rhs.small_list_freq_bytesA[i];
            this->small_list_doc_postingsA[i] += rhs.small_list_doc_postingsA[i];
            this->small_list_freq_postingsA[i] += rhs.small_list_freq_postingsA[i];
        }
        return *this;
    }
};

std::ostream& operator<<(std::ostream& o, const ans_block_size_stats& stats)
{
    o << "full_overhead_bytes = " << stats.full_overhead_bytes << "\n";
    o << "last_overhead_bytes = " << stats.last_overhead_bytes << "\n";
    o << "small_overhead_bytes = " << stats.small_overhead_bytes << "\n";
    o << "total_overhead_bytes = " << stats.total_overhead_bytes << "\n";
    o << "total_postings = " << stats.total_postings << "\n";
    o << "full_block_postings = " << stats.full_block_postings << "\n";
    o << "small_list_postings = " << stats.small_list_postings << "\n";
    o << "last_nonfull_postings = " << stats.last_nonfull_postings << "\n";
    o << "full_block_percent = " << double(stats.full_block_postings) / double(stats.total_postings) * 100 << "\n";
    o << "small_list_percent = " << double(stats.small_list_postings) / double(stats.total_postings) * 100 << "\n";
    o << "last_nonfull_percent = " << double(stats.last_nonfull_postings) / double(stats.total_postings) * 100 << "\n";
    o << "total_doc_bytes = " << stats.total_doc_bytes << "\n";
    o << "total_freq_bytes = " << stats.total_freq_bytes << "\n";
    o << "full_doc_bytes = " << stats.full_doc_bytes << "\n";
    o << "full_freq_bytes = " << stats.full_freq_bytes << "\n";
    o << "small_list_doc_bytes = " << stats.small_list_doc_bytes << "\n";
    o << "small_list_freq_bytes = " << stats.small_list_freq_bytes << "\n";
    o << "last_nonfull_doc_bytes = " << stats.last_nonfull_doc_bytes << "\n";
    o << "last_nonfull_freq_bytes = " << stats.last_nonfull_freq_bytes << "\n";
    o << "full_block_space_percent_docs = " << double(stats.full_doc_bytes) / double(stats.total_doc_bytes) * 100 << "\n";
    o << "full_block_space_percent_freq = " << double(stats.full_freq_bytes) / double(stats.total_freq_bytes) * 100 << "\n";
    o << "small_list_space_percent_docs = " << double(stats.small_list_doc_bytes) / double(stats.total_doc_bytes) * 100 << "\n";
    o << "small_list_space_percent_freq = " << double(stats.small_list_freq_bytes) / double(stats.total_freq_bytes) * 100 << "\n";
    o << "last_nonfull_space_percent_docs = " << double(stats.last_nonfull_doc_bytes) / double(stats.total_doc_bytes) * 100 << "\n";
    o << "last_nonfull_space_percent_freq = " << double(stats.last_nonfull_freq_bytes) / double(stats.total_freq_bytes) * 100 << "\n";
    o << "total_docs_BPI = " << double(stats.total_doc_bytes * 8) / double(stats.total_postings) << "\n";
    o << "total_freqs_BPI = " << double(stats.total_freq_bytes * 8) / double(stats.total_postings) << "\n";
    o << "full_block_docs_BPI = " << double(stats.full_doc_bytes * 8) / double(stats.full_block_postings) << "\n";
    o << "full_block_freqs_BPI = " << double(stats.full_freq_bytes * 8) / double(stats.full_block_postings) << "\n";
    o << "small_list_docs_BPI = " << double(stats.small_list_doc_bytes * 8) / double(stats.small_list_postings) << "\n";
    o << "small_list_freqs_BPI = " << double(stats.small_list_freq_bytes * 8) / double(stats.small_list_postings) << "\n";
    o << "last_nonfull_docs_BPI = " << double(stats.last_nonfull_doc_bytes * 8) / double(stats.last_nonfull_postings) << "\n";
    o << "last_nonfull_freqs_BPI = " << double(stats.last_nonfull_freq_bytes * 8) / double(stats.last_nonfull_postings) << "\n";
    for (size_t i = 1; i < ans::constants::BLOCK_SIZE; i++) {
        if (stats.small_list_doc_postingsA[i]) {
            o << "small_list_docs_BPI[" << i << "] = " << double(stats.small_list_doc_bytesA[i] * 8) / double(stats.small_list_doc_postingsA[i]) << " (";
            o << double(stats.small_list_doc_postingsA[i]) / double(stats.small_list_postings) * 100 << " % small postings - ";
            o << double(stats.small_list_doc_bytesA[i]) / double(stats.small_list_doc_bytes) * 100 << " % small postings space)\n";
        }
    }
    for (size_t i = 1; i < ans::constants::BLOCK_SIZE; i++) {
        if (stats.small_list_freq_postingsA[i]) {
            o << "small_list_freq_BPI[" << i << "] = " << double(stats.small_list_freq_bytesA[i] * 8) / double(stats.small_list_freq_postingsA[i]) << " (";
            o << double(stats.small_list_freq_postingsA[i]) / double(stats.small_list_postings) * 100 << " % small postings - ";
            o << double(stats.small_list_freq_bytesA[i]) / double(stats.small_list_freq_bytes) * 100 << " % small postings space)\n";
        }
    }
    return o << std::dec;
}

template <typename t_ansmodel>
struct ans_block_posting_list {

    template <typename DocsIterator, typename FreqsIterator>
    static void model(std::vector<uint8_t>& doc_model, std::vector<uint8_t>& freq_model, uint32_t n,
        DocsIterator docs_begin, FreqsIterator freqs_begin)
    {
        uint64_t block_size = t_ansmodel::block_size;
        uint64_t blocks = succinct::util::ceil_div(n, block_size);
        DocsIterator docs_it(docs_begin);
        FreqsIterator freqs_it(freqs_begin);
        std::vector<uint32_t> docs_buf(block_size);
        std::vector<uint32_t> freqs_buf(block_size);
        uint32_t last_doc(-1);
        uint32_t block_base = 0;

        for (size_t b = 0; b < blocks; ++b) {
            uint32_t cur_block_size = ((b + 1) * block_size <= n)
                ? block_size
                : (n % block_size);

            for (size_t i = 0; i < cur_block_size; ++i) {
                uint32_t doc(*docs_it++);
                docs_buf[i] = doc - last_doc - 1;
                last_doc = doc;
                freqs_buf[i] = *freqs_it++ - 1;
            }

            t_ansmodel::model(doc_model, docs_buf.data(), last_doc - block_base - (cur_block_size - 1), cur_block_size);
            t_ansmodel::model(freq_model, freqs_buf.data(), uint32_t(-1), cur_block_size);

            block_base = last_doc + 1;
        }
    }

    template <typename DocsIterator, typename FreqsIterator>
    static void write(std::vector<uint8_t>& doc_model, std::vector<uint8_t>& freq_model, std::vector<uint8_t>& out, uint32_t n,
        DocsIterator docs_begin, FreqsIterator freqs_begin)
    {
        TightVariableByte::encode_single(n, out);

        uint64_t block_size = t_ansmodel::block_size;
        uint64_t blocks = succinct::util::ceil_div(n, block_size);
        size_t begin_block_maxs = out.size();
        size_t begin_block_endpoints = begin_block_maxs + 4 * blocks;
        size_t begin_blocks = begin_block_endpoints + 4 * (blocks - 1);
        out.resize(begin_blocks);

        DocsIterator docs_it(docs_begin);
        FreqsIterator freqs_it(freqs_begin);
        std::vector<uint32_t> docs_buf(block_size);
        std::vector<uint32_t> freqs_buf(block_size);
        uint32_t last_doc(-1);
        uint32_t block_base = 0;

        for (size_t b = 0; b < blocks; ++b) {
            uint32_t cur_block_size = ((b + 1) * block_size <= n)
                ? block_size
                : (n % block_size);

            for (size_t i = 0; i < cur_block_size; ++i) {
                uint32_t doc(*docs_it++);
                docs_buf[i] = doc - last_doc - 1;
                last_doc = doc;

                freqs_buf[i] = *freqs_it++ - 1;
            }

            *((uint32_t*)&out[begin_block_maxs + 4 * b]) = last_doc;
            t_ansmodel::encode(docs_buf.data(), last_doc - block_base - (cur_block_size - 1),
                cur_block_size, out, doc_model);

            t_ansmodel::encode(freqs_buf.data(), uint32_t(-1), cur_block_size, out, freq_model);

            if (b != blocks - 1) {
                *((uint32_t*)&out[begin_block_endpoints + 4 * b]) = out.size() - begin_blocks;
            }
            block_base = last_doc + 1;
        }
    }

    class document_enumerator {
    public:
        document_enumerator(uint8_t const* docmodel, uint8_t const* freqmodel, uint8_t const* data, uint64_t universe)
            : m_n(0) // just to silence warnings
            , m_base(TightVariableByte::decode(data, &m_n, 1))
            , m_blocks(succinct::util::ceil_div(m_n, t_ansmodel::block_size))
            , m_block_maxs(m_base)
            , m_block_endpoints(m_block_maxs + 4 * m_blocks)
            , m_blocks_data(m_block_endpoints + 4 * (m_blocks - 1))
            , m_universe(universe)
            , m_doc_model_data(docmodel)
            , m_freqs_model_data(freqmodel)
        {
            m_docs_buf.resize(t_ansmodel::block_size);
            m_freqs_buf.resize(t_ansmodel::block_size);
            reset();
        }

        void reset()
        {
            decode_docs_block(0);
        }

        void QS_ALWAYSINLINE next()
        {
            ++m_pos_in_block;
            if (QS_UNLIKELY(m_pos_in_block == m_cur_block_size)) {
                if (m_cur_block + 1 == m_blocks) {
                    m_cur_docid = m_universe;
                    return;
                }
                decode_docs_block(m_cur_block + 1);
            } else {
                m_cur_docid += m_docs_buf[m_pos_in_block] + 1;
            }
        }

        void QS_ALWAYSINLINE next_geq(uint64_t lower_bound)
        {
            assert(lower_bound >= m_cur_docid);
            if (QS_UNLIKELY(lower_bound > m_cur_block_max)) {
                // binary search seems to perform worse here
                if (lower_bound > block_max(m_blocks - 1)) {
                    m_cur_docid = m_universe;
                    return;
                }

                uint64_t block = m_cur_block + 1;
                while (block_max(block) < lower_bound) {
                    ++block;
                }

                decode_docs_block(block);
            }

            while (docid() < lower_bound) {
                m_cur_docid += m_docs_buf[++m_pos_in_block] + 1;
                assert(m_pos_in_block < m_cur_block_size);
            }
        }

        void QS_ALWAYSINLINE move(uint64_t pos)
        {
            assert(pos >= position());
            uint64_t block = pos / t_ansmodel::block_size;
            if (QS_UNLIKELY(block != m_cur_block)) {
                decode_docs_block(block);
            }
            while (position() < pos) {
                m_cur_docid += m_docs_buf[++m_pos_in_block] + 1;
            }
        }

        uint64_t docid() const
        {
            return m_cur_docid;
        }

        uint64_t QS_ALWAYSINLINE freq()
        {
            if (!m_freqs_decoded) {
                decode_freqs_block();
            }
            return m_freqs_buf[m_pos_in_block] + 1;
        }

        uint64_t position() const
        {
            return m_cur_block * t_ansmodel::block_size + m_pos_in_block;
        }

        uint64_t size() const
        {
            return m_n;
        }

        ans_block_size_stats size_stats() const
        {
            ans_block_size_stats bss;

            uint8_t const* doc_ptr = m_blocks_data;
            static const uint64_t block_size = t_ansmodel::block_size;
            std::vector<uint32_t> buf(block_size);
            for (size_t b = 0; b < m_blocks; ++b) {
                uint32_t cur_block_size = ((b + 1) * block_size <= size())
                    ? block_size
                    : (size() % block_size);

                uint32_t cur_base = (b ? block_max(b - 1) : uint32_t(-1)) + 1;
                uint8_t const* freq_ptr = t_ansmodel::decode(doc_ptr, buf.data(),
                    block_max(b) - cur_base - (cur_block_size - 1),
                    cur_block_size, m_doc_model_data);
                uint8_t const* end_ptr = t_ansmodel::decode(freq_ptr, buf.data(),
                    uint32_t(-1), cur_block_size, m_freqs_model_data);

                size_t doc_bytes = freq_ptr - doc_ptr;
                size_t freq_bytes = end_ptr - freq_ptr;
                doc_ptr = end_ptr;

                if (b == 0) {
                    bss.total_overhead_bytes += sizeof(uint32_t); // block max only
                } else {
                    bss.total_overhead_bytes += 2 * sizeof(uint32_t); // block max and data offset
                }

                bss.total_doc_bytes += doc_bytes;
                bss.total_freq_bytes += freq_bytes;
                bss.total_postings += cur_block_size;
                if (m_n < block_size) {
                    bss.small_overhead_bytes += sizeof(uint32_t); // block max only
                    bss.small_list_doc_bytes = doc_bytes;
                    bss.small_list_freq_bytes = freq_bytes;
                    bss.small_list_postings = m_n;
                    bss.small_list_doc_postingsA[m_n] += m_n;
                    bss.small_list_freq_postingsA[m_n] += m_n;
                    bss.small_list_doc_bytesA[m_n] += doc_bytes;
                    bss.small_list_freq_bytesA[m_n] += freq_bytes;
                } else {
                    if (cur_block_size < block_size) {
                        bss.last_overhead_bytes += 2 * sizeof(uint32_t); // block max and data offset
                        bss.last_nonfull_doc_bytes = doc_bytes;
                        bss.last_nonfull_freq_bytes = freq_bytes;
                        bss.last_nonfull_postings = cur_block_size;
                    } else {
                        if (b == 0) {
                            bss.full_overhead_bytes += sizeof(uint32_t); // block max and data offset
                        } else {
                            bss.full_overhead_bytes += 2 * sizeof(uint32_t); // block max and data offset
                        }
                        bss.full_doc_bytes += doc_bytes;
                        bss.full_freq_bytes += freq_bytes;
                        bss.full_block_postings += block_size;
                    }
                }
            }

            return bss;
        }

        uint64_t stats_freqs_size() const
        {
            uint64_t bytes = 0;
            uint8_t const* ptr = m_blocks_data;
            static const uint64_t block_size = t_ansmodel::block_size;
            std::vector<uint32_t> buf(block_size);
            for (size_t b = 0; b < m_blocks; ++b) {
                uint32_t cur_block_size = ((b + 1) * block_size <= size())
                    ? block_size
                    : (size() % block_size);

                uint32_t cur_base = (b ? block_max(b - 1) : uint32_t(-1)) + 1;
                uint8_t const* freq_ptr = t_ansmodel::decode(ptr, buf.data(),
                    block_max(b) - cur_base - (cur_block_size - 1),
                    cur_block_size, m_doc_model_data);
                ptr = t_ansmodel::decode(freq_ptr, buf.data(),
                    uint32_t(-1), cur_block_size, m_freqs_model_data);
                bytes += ptr - freq_ptr;
            }

            return bytes;
        }

    private:
        uint32_t block_max(uint32_t block) const
        {
            return ((uint32_t const*)m_block_maxs)[block];
        }

        void QS_NOINLINE decode_docs_block(uint64_t block)
        {

            static const uint64_t block_size = t_ansmodel::block_size;
            uint32_t endpoint = block
                ? ((uint32_t const*)m_block_endpoints)[block - 1]
                : 0;
            uint8_t const* block_data = m_blocks_data + endpoint;

            m_cur_block_size = ((block + 1) * block_size <= size())
                ? block_size
                : (size() % block_size);
            uint32_t cur_base = (block ? block_max(block - 1) : uint32_t(-1)) + 1;
            m_cur_block_max = block_max(block);
            m_freqs_block_data = t_ansmodel::decode(block_data, m_docs_buf.data(),
                m_cur_block_max - cur_base - (m_cur_block_size - 1),
                m_cur_block_size, m_doc_model_data);

            m_docs_buf[0] += cur_base;

            m_cur_block = block;
            m_pos_in_block = 0;
            m_cur_docid = m_docs_buf[0];
            m_freqs_decoded = false;
        }

        void QS_NOINLINE decode_freqs_block()
        {
            t_ansmodel::decode(m_freqs_block_data, m_freqs_buf.data(),
                uint32_t(-1), m_cur_block_size, m_freqs_model_data);
            m_freqs_decoded = true;
        }

        uint32_t m_n;
        uint8_t const* m_base;
        uint32_t m_blocks;
        uint8_t const* m_block_maxs;
        uint8_t const* m_block_endpoints;
        uint8_t const* m_blocks_data;
        uint64_t m_universe;

        uint32_t m_cur_block;
        uint32_t m_pos_in_block;
        uint32_t m_cur_block_max;
        uint32_t m_cur_block_size;
        uint32_t m_cur_docid;

        uint8_t const* m_freqs_block_data;
        bool m_freqs_decoded;

        std::vector<uint32_t> m_docs_buf;
        std::vector<uint32_t> m_freqs_buf;

        uint8_t const* m_doc_model_data;
        uint8_t const* m_freqs_model_data;
    };
};
}
