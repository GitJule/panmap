#include <random>
#include <seqan3/io/sequence_file/output.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <filesystem>
#include <loadjst.hpp>

void generate_test_sequences(const rcs_store_t& jst_data, 
                            const std::filesystem::path& output_path,
                            size_t pattern_length = 50,
                            size_t samples_per_node = 1000,
                            unsigned seed = 42) 
{
    std::mt19937 gen(seed);
    seqan3::sequence_file_output fasta_out{output_path};

    auto search_tree = libjst::make_volatile(jst_data)
        | libjst::labelled()
        | libjst::trim(pattern_length - 1)
        | libjst::merge(libjst::merge_adjacent)
        | libjst::seek();

    libjst::tree_traverser_base traverser{search_tree};

    for (auto const& node : traverser) {
        auto const& seq = node.sequence();
        if (seq.empty() || seq.size() < pattern_length)
            continue;

        std::uniform_int_distribution<size_t> distrib(0, seq.size() - pattern_length);

        for (size_t i = 0; i < samples_per_node; ++i) {
            size_t start = distrib(gen);
            auto safe_slice = seq | seqan3::views::slice(start, start + pattern_length);

            fasta_out.emplace_back(
                safe_slice,
                std::to_string(node.position().variant_id()) + "_" + std::to_string(start)
            );
        }
    }
}
