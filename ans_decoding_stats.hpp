#pragma once

#include <cstdint>
#include <iostream>

struct ans_dec_stats {
    uint64_t total_postings_decodes = 0;
    uint64_t num_no_decodes = 0;
    uint64_t postings_no_decodes = 0;
    uint64_t num_model_decodes = 0;
    uint64_t postings_model_decodes = 0;
    uint64_t model_usage[256] = { 0 };
    uint64_t model_frame_size[256] = { 0 };
    uint64_t final_state_bytes[8] = { 0 };
    ans_dec_stats()
    {
    }
    static ans_dec_stats& doc_stats()
    {
        static std::unique_ptr<ans_dec_stats> s(new ans_dec_stats);
        return *s.get();
    }
    static ans_dec_stats& freq_stats()
    {
        static std::unique_ptr<ans_dec_stats> s(new ans_dec_stats);
        return *s.get();
    }
    static ans_dec_stats& stats(bool is_freq)
    {
        if (is_freq)
            return freq_stats();
        return doc_stats();
    }
    static void print_doc_stats()
    {
        std::cout << "DOCUMENT DECODING STATS" << std::endl;
        const auto& s = doc_stats();
        std::cout << "total_postings_decoded = " << s.total_postings_decodes << std::endl;
        std::cout << "num_no_decodes = " << s.num_no_decodes << std::endl;
        std::cout << "postings_no_decodes = " << s.postings_no_decodes << std::endl;
        std::cout << "num_model_decodes = " << s.num_model_decodes << std::endl;
        std::cout << "postings_model_decodes = " << s.postings_model_decodes << std::endl;
        size_t non_zero = 0;
        double avg_model = 0;
        double total_model = 0;
        for (size_t i = 0; i < 256; i++) {
            if (s.model_usage[i] != 0) {
                std::cout << non_zero++ << " model_id = " << i << " M = " << s.model_frame_size[i]
                          << " usage = " << s.model_usage[i]
                          << " percent = " << 100.0 * double(s.model_usage[i]) / double(s.num_model_decodes)
                          << std::endl;
                avg_model += ((i + 1) * s.model_usage[i]);
                total_model += s.model_usage[i];
            }
        }
        std::cout << "average model_id = " << (avg_model / total_model) - 1.0 << std::endl;
        for (size_t i = 0; i < 8; i++) {
            std::cout << "final_state_bytes[" << i << "] = " << s.final_state_bytes[i]
                      << " percent = " << 100.0 * double(s.final_state_bytes[i]) / double(s.num_model_decodes)
                      << std::endl;
        }
    }
    static void print_freq_stats()
    {
        std::cout << "FREQ DECODING STATS" << std::endl;
        const auto& s = freq_stats();
        std::cout << "total_postings_decoded = " << s.total_postings_decodes << std::endl;
        std::cout << "num_no_decodes = " << s.num_no_decodes << std::endl;
        std::cout << "postings_no_decodes = " << s.postings_no_decodes << std::endl;
        std::cout << "num_model_decodes = " << s.num_model_decodes << std::endl;
        std::cout << "postings_model_decodes = " << s.postings_model_decodes << std::endl;
        size_t non_zero = 0;
        double avg_model = 0;
        double total_model = 0;
        for (size_t i = 0; i < 256; i++) {
            if (s.model_usage[i] != 0) {
                std::cout << non_zero++ << "model_id = " << i << " M = " << s.model_frame_size[i]
                          << " usage = " << s.model_usage[i]
                          << " percent = " << 100.0 * double(s.model_usage[i]) / double(s.num_model_decodes)
                          << std::endl;
                avg_model += ((i + 1) * s.model_usage[i]);
                total_model += s.model_usage[i];
            }
        }
        std::cout << "average model_id = " << (avg_model / total_model) - 1.0 << std::endl;
        for (size_t i = 0; i < 8; i++) {
            std::cout << "final_state_bytes[" << i << "] = " << s.final_state_bytes[i]
                      << " percent = " << 100.0 * double(s.final_state_bytes[i]) / double(s.num_model_decodes)
                      << std::endl;
        }
    }
};
