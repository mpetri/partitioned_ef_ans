#pragma once

#include <succinct/bit_vector.hpp>
#include <succinct/mappable_vector.hpp>

#include "ans_block_posting_list.hpp"
#include "compact_elias_fano.hpp"

#include "ans_models.hpp"

namespace quasi_succinct {

template <typename t_ansmodel>
class ans_block_freq_index {
public:
    ans_block_freq_index()
        : m_size(0)
    {
    }

    class builder {
    public:
        builder(uint64_t num_docs, global_parameters const& params)
            : m_params(params)
        {
            m_num_docs = num_docs;
            m_endpoints.push_back(0);
            m_doc_counts = t_ansmodel::create_empty_counts();
            m_freq_counts = t_ansmodel::create_empty_counts();
        }

        template <typename DocsIterator, typename FreqsIterator>
        void model_posting_list(uint64_t n, DocsIterator docs_begin,
            FreqsIterator freqs_begin, uint64_t /*occurrences*/)
        {
            ans_block_posting_list<t_ansmodel>::model(m_doc_counts, m_freq_counts, n, docs_begin, freqs_begin);
        }

        template <typename DocsIterator, typename FreqsIterator>
        void add_posting_list(uint64_t n, DocsIterator docs_begin,
            FreqsIterator freqs_begin, uint64_t /* occurrences */)
        {
            if (!n)
                throw std::invalid_argument("List must be nonempty");
            ans_block_posting_list<t_ansmodel>::write(m_doc_enc_model, m_freq_enc_model, m_lists, n, docs_begin, freqs_begin);
            m_endpoints.push_back(m_lists.size());
        }

        void freeze_models()
        {
            m_doc_enc_model = t_ansmodel::create_enc_model_from_counts(m_doc_counts);
            m_freq_enc_model = t_ansmodel::create_enc_model_from_counts(m_freq_counts);
        }

        void build(ans_block_freq_index& sq)
        {
            sq.m_params = m_params;
            sq.m_size = m_endpoints.size() - 1;
            sq.m_num_docs = m_num_docs;
            sq.m_lists.steal(m_lists);
            sq.m_doc_dec_model = t_ansmodel::create_dec_model(m_doc_enc_model);
            sq.m_freq_dec_model = t_ansmodel::create_dec_model(m_freq_enc_model);
            m_doc_enc_model.clear();
            m_freq_enc_model.clear();

            succinct::bit_vector_builder bvb;
            compact_elias_fano::write(bvb, m_endpoints.begin(),
                sq.m_lists.size(), sq.m_size,
                m_params); // XXX
            succinct::bit_vector(&bvb).swap(sq.m_endpoints);
        }

    private:
        global_parameters m_params;
        size_t m_num_docs;
        std::vector<uint64_t> m_endpoints;
        std::vector<uint8_t> m_lists;
        std::vector<uint8_t> m_doc_counts;
        std::vector<uint8_t> m_freq_counts;
        std::vector<uint8_t> m_doc_enc_model;
        std::vector<uint8_t> m_freq_enc_model;
        std::vector<uint8_t> m_doc_dec_model;
        std::vector<uint8_t> m_freq_dec_model;
    };

    size_t size() const
    {
        return m_size;
    }

    uint64_t num_docs() const
    {
        return m_num_docs;
    }

    typedef typename ans_block_posting_list<t_ansmodel>::document_enumerator document_enumerator;

    document_enumerator operator[](size_t i) const
    {
        assert(i < size());
        compact_elias_fano::enumerator endpoints(m_endpoints, 0,
            m_lists.size(), m_size,
            m_params);

        auto endpoint = endpoints.move(i).second;
        return document_enumerator(m_doc_dec_model.data(), m_freq_dec_model.data(), m_lists.data() + endpoint, num_docs());
    }

    void swap(ans_block_freq_index& other)
    {
        std::swap(m_params, other.m_params);
        std::swap(m_size, other.m_size);
        m_endpoints.swap(other.m_endpoints);
        m_lists.swap(other.m_lists);
        m_doc_dec_model.swap(other.m_doc_dec_model);
        m_freq_dec_model.swap(other.m_freq_dec_model);
    }

    template <typename Visitor>
    void map(Visitor& visit)
    {
        // clang-format off
        visit(m_params, "m_params")
            (m_size, "m_size")
            (m_num_docs, "m_num_docs")
            (m_endpoints, "m_endpoints")
            (m_lists, "m_lists")
            (m_doc_dec_model, "m_doc_dec_model")
            (m_freq_dec_model, "m_freq_dec_model");
        // clang-format on
    }

private:
    global_parameters m_params;
    size_t m_size;
    size_t m_num_docs;
    succinct::bit_vector m_endpoints;
    succinct::mapper::mappable_vector<uint8_t> m_lists;
    succinct::mapper::mappable_vector<uint8_t> m_doc_dec_model;
    succinct::mapper::mappable_vector<uint8_t> m_freq_dec_model;
};
}
