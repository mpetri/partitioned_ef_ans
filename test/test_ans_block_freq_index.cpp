#define BOOST_TEST_MODULE ans_block_freq_index

#include "test_generic_sequence.hpp"

#include "ans_block_freq_index.hpp"
#include "ans_packed_model.hpp"
#include <succinct/mapper.hpp>

#include <algorithm>
#include <cstdlib>
#include <random>
#include <vector>

template <typename t_ans_model>
void test_ans_block_freq_index_small()
{
    quasi_succinct::global_parameters params;
    uint64_t universe = 20000;
    typedef quasi_succinct::ans_block_freq_index<t_ans_model> collection_type;
    typename collection_type::builder b(universe, params);

    typedef std::vector<uint64_t> vec_type;
    std::vector<std::pair<vec_type, vec_type>> posting_lists(300);
    for (auto& plist : posting_lists) {
        double avg_gap = 1.1 + double(rand()) / RAND_MAX * 10;
        uint64_t n = uint64_t(universe / avg_gap);
        plist.first = random_sequence(universe, n, true);
        plist.second.resize(n);
        std::generate(plist.second.begin(), plist.second.end(),
            []() { return (rand() % 256) + 1; });

        b.model_posting_list(n, plist.first.begin(),
            plist.second.begin(), 0);
    }
    b.freeze_models();

    for (auto& plist : posting_lists) {
        b.add_posting_list(plist.second.size(), plist.first.begin(),
            plist.second.begin(), 0);
    }

    {
        collection_type coll;
        b.build(coll);
        succinct::mapper::freeze(coll, "temp.bin");
    }

    {
        collection_type coll;
        boost::iostreams::mapped_file_source m("temp.bin");
        succinct::mapper::map(coll, m);

        for (size_t i = 0; i < posting_lists.size(); ++i) {
            auto const& plist = posting_lists[i];
            auto doc_enum = coll[i];
            BOOST_REQUIRE_EQUAL(plist.first.size(), doc_enum.size());
            for (size_t p = 0; p < plist.first.size(); ++p, doc_enum.next()) {
                MY_REQUIRE_EQUAL(plist.first[p], doc_enum.docid(),
                    "i = " << i << " p = " << p << "/" << plist.first.size());
                MY_REQUIRE_EQUAL(plist.second[p], doc_enum.freq(),
                    "i = " << i << " n = " << plist.first.size() << " p = " << p << " p%128 = " << p % 128);
            }
            BOOST_REQUIRE_EQUAL(coll.num_docs(), doc_enum.docid());
        }
    }
}

template <typename t_ans_model>
void test_ans_block_freq_index_large(size_t avg_gap = 10)
{
    quasi_succinct::global_parameters params;
    uint64_t universe = 26000000;
    typedef quasi_succinct::ans_block_freq_index<t_ans_model> collection_type;
    typename collection_type::builder b(universe, params);

    typedef std::vector<uint64_t> vec_type;
    std::vector<std::pair<vec_type, vec_type>> posting_lists(20);

    std::mt19937 gen(rand());
    std::normal_distribution<> dis(avg_gap, 2);
    for (auto& plist : posting_lists) {
        double avg_gap = 1.1 + dis(gen);
        uint64_t n = uint64_t(universe / avg_gap);
        plist.first = random_sequence(universe, n, true);
        plist.second.resize(n);
        std::generate(plist.second.begin(), plist.second.end(),
            []() { return (rand() % 256) + 1; });

        b.model_posting_list(n, plist.first.begin(),
            plist.second.begin(), 0);
    }
    b.freeze_models();

    for (auto& plist : posting_lists) {
        b.add_posting_list(plist.second.size(), plist.first.begin(),
            plist.second.begin(), 0);
    }

    {
        collection_type coll;
        b.build(coll);
        succinct::mapper::freeze(coll, "temp.bin");
    }

    {
        collection_type coll;
        boost::iostreams::mapped_file_source m("temp.bin");
        succinct::mapper::map(coll, m);

        for (size_t i = 0; i < posting_lists.size(); ++i) {
            auto const& plist = posting_lists[i];
            auto doc_enum = coll[i];
            BOOST_REQUIRE_EQUAL(plist.first.size(), doc_enum.size());
            for (size_t p = 0; p < plist.first.size(); ++p, doc_enum.next()) {
                MY_REQUIRE_EQUAL(plist.first[p], doc_enum.docid(),
                    "i = " << i << " n = " << p << " / " << plist.first.size() << p << " p%128 = " << p % 128);
                MY_REQUIRE_EQUAL(plist.second[p], doc_enum.freq(),
                    "i = " << i << " n = " << plist.first.size() << " p = " << p << " p%128 = " << p % 128);
            }
            BOOST_REQUIRE_EQUAL(coll.num_docs(), doc_enum.docid());
        }
    }
}

BOOST_AUTO_TEST_CASE(ans_msb_model_model_max_1d_small)
{
    test_ans_block_freq_index_small<quasi_succinct::ans_msb_model<msb_model_max_1d>>();
}

BOOST_AUTO_TEST_CASE(ans_msb_model_model_max_1d_large)
{
    // test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_max_1d>>(10);
    // test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_max_1d>>(100);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_max_1d>>(1000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_max_1d>>(100000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_max_1d>>(1000000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_max_1d>>(10000000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_max_1d>>(20000000);
}

BOOST_AUTO_TEST_CASE(ans_msb_model_model_minmax_2d_small)
{
    test_ans_block_freq_index_small<quasi_succinct::ans_msb_model<msb_model_minmax_2d>>();
}

BOOST_AUTO_TEST_CASE(ans_msb_model_model_minmax_2d_large)
{
    // test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_minmax_2d>>(10);
    // test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_minmax_2d>>(100);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_minmax_2d>>(1000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_minmax_2d>>(100000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_minmax_2d>>(1000000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_minmax_2d>>(10000000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_minmax_2d>>(20000000);
}

BOOST_AUTO_TEST_CASE(ans_msb_model_model_med90p_2d_large)
{
    // test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d>>(10);
    // test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d>>(100);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d>>(1000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d>>(100000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d>>(1000000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d>>(10000000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d>>(20000000);
}

BOOST_AUTO_TEST_CASE(ans_msb_model_model_med90p_es_2d_large)
{
    // test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d_es>>(10);
    // test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d_es>>(100);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d_es>>(1000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d_es>>(100000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d_es>>(1000000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d_es>>(10000000);
    test_ans_block_freq_index_large<quasi_succinct::ans_msb_model<msb_model_med90p_2d_es>>(20000000);
}

// BOOST_AUTO_TEST_CASE(ans_packed_model_model_max_1d_small)
// {
//     test_ans_block_freq_index_small<quasi_succinct::ans_packed_model<model_max_1d>>();
// }

// BOOST_AUTO_TEST_CASE(ans_packed_model_model_minmax_2d_small)
// {
//     test_ans_block_freq_index_small<quasi_succinct::ans_packed_model<model_minmax_2d>>();
// }

// BOOST_AUTO_TEST_CASE(ans_packed_model_model_med90p_2d_small)
// {
//     test_ans_block_freq_index_small<quasi_succinct::ans_packed_model<model_med90p_2d>>();
// }

// BOOST_AUTO_TEST_CASE(ans_packed_model_model_max_1d_large)
// {
//     test_ans_block_freq_index_large<quasi_succinct::ans_packed_model<model_max_1d>>();
// }

// BOOST_AUTO_TEST_CASE(ans_packed_model_model_minmax_2d_large)
// {
//     test_ans_block_freq_index_large<quasi_succinct::ans_packed_model<model_minmax_2d>>();
// }

// BOOST_AUTO_TEST_CASE(ans_packed_model_model_med90p_2d_large)
// {
//     test_ans_block_freq_index_large<quasi_succinct::ans_packed_model<model_med90p_2d>>();
// }
