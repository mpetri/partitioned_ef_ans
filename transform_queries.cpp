
#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>

#include "libstemmer.h"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

po::variables_map parse_cmdargs(int argc, char const* argv[])
{
    po::variables_map vm;
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help,h", "produce help message")
        ("query-file,q",po::value<std::string>()->required(), "input query file")
        ("term-file,t",po::value<std::string>()->required(), "collection term->id file");
    // clang-format on
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            std::cout << desc << "\n";
            exit(EXIT_SUCCESS);
        }
        po::notify(vm);
    } catch (const po::required_option& e) {
        std::cerr << desc;
        std::cerr << "Missing required option: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    } catch (po::error& e) {
        std::cerr << desc;
        std::cerr << "Error parsing cmdargs: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    return vm;
}

std::vector<std::vector<std::string>>
parse_query_file(std::string query_file)
{
    std::vector<std::vector<std::string>> unstemmed_queries;
    std::ifstream qfs(query_file);
    if (qfs) {
        for (std::string query; std::getline(qfs, query);) {
            std::istringstream q(query);
            std::vector<std::string> tokenized_qry;
            for (std::string tok; std::getline(q, tok, ' ');) {
                if (tok.size())
                    tokenized_qry.push_back(tok);
            }
            if (tokenized_qry.size()) {
                unstemmed_queries.push_back(tokenized_qry);
            }
        }
    }
    std::cerr << "Parsed " << unstemmed_queries.size() << " queries" << std::endl;
    return unstemmed_queries;
}

std::unordered_map<std::string, uint64_t>
parse_term_file(std::string term_file)
{
    std::unordered_map<std::string, uint64_t> tm;
    std::ifstream tfs(term_file);
    if (tfs) {
        uint64_t cur_id = 0;
        for (std::string term; std::getline(tfs, term);) {
            tm[term] = cur_id++;
        }
    }
    std::cerr << "Loaded " << tm.size() << " terms" << std::endl;
    return tm;
}

std::vector<std::vector<std::string>>
stem_queries(const std::vector<std::vector<std::string>>& unstemmed_queries)
{
    struct sb_stemmer* stemmer = sb_stemmer_new("porter", NULL);
    if (stemmer == 0) {
        std::cerr << "Error creating stemmer" << std::endl;
        exit(1);
    }

    std::vector<std::vector<std::string>> stemmed_queries;
    for (const auto& qry : unstemmed_queries) {
        std::vector<std::string> sq;
        for (const auto& tok : qry) {
            std::string ttok = tok;
            std::transform(ttok.begin(), ttok.end(), ttok.begin(), ::tolower);
            const char* stemmed = (const char*)sb_stemmer_stem(stemmer, (const sb_symbol*)ttok.c_str(), ttok.size());
            std::string stok(stemmed);
            if (!stok.empty()) {
                sq.push_back(stok);
            }
        }
        if (sq.size()) {
            stemmed_queries.push_back(sq);
        }
    }

    sb_stemmer_delete(stemmer);
    return stemmed_queries;
}

std::vector<std::vector<uint64_t>>
map_queries(const std::vector<std::vector<std::string>>& stemmed_queries,
    const std::unordered_map<std::string, uint64_t>& tm)
{
    uint64_t num_skipped = 0;
    std::vector<std::vector<uint64_t>> mapped_queries;
    for (const auto& qry : stemmed_queries) {
        std::vector<uint64_t> mq;
        bool skip = false;
        for (const auto& tok : qry) {
            auto itr = tm.find(tok);
            if (itr == tm.end()) {
                std::cerr << "Could not map token '" << tok << "'. skipping query." << std::endl;
                skip = true;
            } else {
                mq.push_back(itr->second);
            }
        }
        if (!skip) {
            mapped_queries.push_back(mq);
        } else {
            num_skipped++;
        }
    }
    std::cerr << "Mapped " << mapped_queries.size() << " to collection ids." << std::endl;
    std::cerr << "Skipped " << num_skipped << " queries due to unknown tokens." << std::endl;
    return mapped_queries;
}

int main(int argc, const char** argv)
{
    auto cmdargs = parse_cmdargs(argc, argv);
    auto query_file = cmdargs["query-file"].as<std::string>();
    auto term_file = cmdargs["term-file"].as<std::string>();

    //(1) parse inputs
    auto unstemmed_queries = parse_query_file(query_file);
    auto term_map = parse_term_file(term_file);

    // (2) stem using porter stemmer
    auto stemmed_queries = stem_queries(unstemmed_queries);

    // (3) map to collection ids
    auto mapped_queries = map_queries(stemmed_queries, term_map);

    // (4) output to stdout
    for (const auto& q : mapped_queries) {
        for (size_t i = 0; i < q.size() - 1; i++)
            std::cout << q[i] << "\t";
        std::cout << q.back() << std::endl;
    }
}
