// // SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// // SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// // SPDX-License-Identifier: CC0-1.0

#include <sharg/all.hpp>
#include <sharg/parser.hpp>
#include <fstream>
#include <span>
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <ranges> 
#include <filesystem>  

#include <configuration.hpp>
#include <loadjst.hpp>
#include <seqan3/search/search.hpp>
#include <libspm/matcher/horspool_matcher.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>
#include <seqan3/utility/views/slice.hpp>

struct match_position
{
    libjst::seek_position jst_seek_position; 
    std::ptrdiff_t label_offset;
};

std::vector<match_position> jst_search(const rcs_store_t& jst_data, const reference_t & read)
{
    spm::horspool_matcher matcher{read};
    size_t const needle_size = spm::window_size(matcher);
    
    auto search_tree = libjst::make_volatile(jst_data)
                 | libjst::labelled()
                 | libjst::trim(needle_size -1)
                 | libjst::left_extend(needle_size -1)
                 | libjst::merge() // make big nodes
                 | libjst::seek();

    std::vector<match_position> hits;
    hits.reserve(1000);
    libjst::tree_traverser_base traverser{search_tree};

    for (auto it = traverser.begin(); it != traverser.end(); ++it) 
    {
        auto && cargo = *it;
        matcher(cargo.sequence(), [&] ([[maybe_unused]] auto && label_finder) {
                hits.emplace_back(
                cargo.position(), 
                static_cast<std::ptrdiff_t>(label_finder.offset()));
            });
    }
    return hits;
}

#include <random>

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
        | libjst::merge()
        | libjst::seek();

    libjst::tree_traverser_base traverser{search_tree};
    
    for (auto const& node : traverser) { 
        auto const& seq = node.sequence();
        if (seq.empty() || seq.size() < pattern_length) continue;

        std::uniform_int_distribution<size_t> distrib(0, seq.size() - pattern_length);
        
        for (size_t i = 0; i < samples_per_node; ++i) {
            size_t start = distrib(gen);
            auto safe_slice = seq | seqan3::views::slice(start, std::min(start + pattern_length, seq.size()));
            
            fasta_out.emplace_back(
                safe_slice,
                std::to_string(jst_data.reference.at(ref_id)) + "_" +
                std::to_string(start)
            );
        }
    }
}

void map_reads(std::filesystem::path const & query_path,
               std::filesystem::path const & sam_path,
               const rcs_store_t & jst_data,
               [[maybe_unused]] uint8_t const errors)
{
    seqan3::sequence_file_input query_file_in{query_path};
    seqan3::sam_file_output sam_out{sam_path}; 

    sam_out.header().ref_dict = jst_data.reference;
    sam_out.header().ref_id_info = jst_data.reference_lengths;

    auto const align_config = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                               seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                                               seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                               seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}} | 
                                                               seqan3::align_cfg::edit_scheme | 
                                                               seqan3::align_cfg::output_alignment{} |
                                                               seqan3::align_cfg::output_begin_position{} |
                                                               seqan3::align_cfg::output_end_position{} |
                                                               seqan3::align_cfg::output_score{};


    for (auto && record : query_file_in) // beginnt Schleife, die jeden Eintrag in der Query-Datei durchläuft
    {
        const auto & query_seq = record.sequence()
        {
            if (query_seq.empty()) continue;
        }
        auto hits = jst_search(jst_data, query_seq_size);
        size_t const window_size = std::ranges::size(query_seq);
        if (window_size == 0) continue;

        for (auto & hit_pos : hits)
        {auto search_tree = libjst::make_volatile(jst_data)
                          | libjst::labelled()
                          | libjst::trim(window_size)
                          | libjst::left_extend(window_size)
                          | libjst::merge() // make big nodes
                          | libjst::seek();

        auto node_covering_hit = search_tree.seek(hit_pos.jst_seek_position);
        auto cargo_hit = *node_covering_hit;
        auto jst_branch_sequence = cargo_hit.sequence();

        std::ptrdiff_t begin_pos = std::max<std::ptrdiff_t>(0, hit_pos.label_offset - errors);
        std::ptrdiff_t end_pos = std::min<std::ptrdiff_t>(hit_pos.label_offset + std::ranges::size(query_seq_size) + errors,
                                                               std::ranges::size(jst_branch_sequence));

        auto reference_segment = jst_branch_sequence | seqan3::views::slice(begin_pos, end_pos);

            
            for (auto&& alignment : seqan3::align_pairwise(std::tie(reference_segment,query_seq), align_config))
            {
                auto cigar = seqan3::cigar_from_alignment(alignment.alignment(), alignment.sequence1_begin_position(), // Startpos Referenz
                alignment.sequence2_begin_position());  // Startpos Query // wandelt das Alignment in einen CIGAR-String um

                auto ref_pos = hit_pos.jst_seek_position.get_variant_position() + begin_pos;

                sam_out.emplace_back(query_seq_size,
                                     record.id(),
                                     jst_data.ids[hit_pos.jst_seek_position],
                                     ref_pos + 1,                                   
                                     cigar,
                                     record.base_qualities(),
                                     60u + alignment.score());
                                     seqan3::sam_tag_dictionary{},
                                     record.sequence(); 
            }
        }
    }
}
 
void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & query_path,
                 std::filesystem::path const & sam_path,
                 uint8_t const errors)
{
 auto jst_data = loadjst(reference_path);
 map_reads(query_path,
           sam_path, 
           jst_data,
           errors); 
}

void initialise_argument_parser(sharg::parser & parser, configuration & args)
{
    parser.info.author = "Julia Bodnar"; // setzt den Autor des Programms
    parser.info.short_description = "JST-based read mapper"; // setzt Beschreibung des Programms 
    parser.info.version = "1.0.0"; // setzt Version des Programms 

    parser.add_option(args.reference_path, // fügt eine Option für den Pfad zur Referenzdatei hinzu
                      sharg::config{.short_id = 'r',
                                    .long_id= "reference", 
                                    .description = "referenzgenome (JST-format)",
                                    .required = true, 
                                    .validator = sharg::input_file_validator{{"jst"}}});
                    
    parser.add_option(args.query_path, 
                      sharg::config{.short_id = 'q',
                                    .long_id= "query", 
                                    .description = "Input reads", 
                                    .required = true, 
                                    .validator = sharg::input_file_validator{{"fq","fastq"}}});

    parser.add_option(args.sam_path, // fügt eine Option für den Pfad zur SAM Ausgabedatei hinzu
                      sharg::config{.short_id = 'o',
                                    .long_id= "output", 
                                    .description = "The output SAM file", 
                                    .required = true, 
                                    .validator = sharg::output_file_validator{{"sam"}}});

    parser.add_option(args.errors, // fügt eine Option für die maximale Anzahl von Feldern hinzu
                      sharg::config{.short_id = 'e',
                                    .long_id= "error", 
                                    .description = "Maximum allowed errors.", 
                                    .required = true, 
                                    .validator = sharg::arithmetic_range_validator{0, 4}});                         
}

int main(int argc, char const** argv)
{
    sharg::parser parser{"PanMap", argc, argv};
    configuration args{};
    
    initialise_argument_parser(parser, args);
    
    run_program(args.reference_path, 
               args.query_path, 
               args.sam_path, 
               args.errors);
    
    return 0;
}
