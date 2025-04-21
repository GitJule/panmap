// // SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// // SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// // SPDX-License-Identifier: CC0-1.0

#include <sharg/all.hpp>
#include <fstream>
#include <span>
#include <ranges>
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>

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


// struct jst_adapter 
// {
//     using iterator = typename std::vector<seqan3::dna5_vector>::const_iterator;
    
//     const JST_Data& data;
    
//     size_t size() const { return data.sequences.size(); }
//     auto operator[](size_t i) const { return std::views::all(data.sequences[i]); }

//     iterator begin() const { return data.sequences.begin(); }
//     iterator end() const { return data.sequences.end(); }
// };

std::vector<size_t> jst_search(const auto& jst_data, const reference_t & read)
{
    // auto config = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}} //keine fehler erlaubt
    //             | seqan3::search_cfg::hit_all_best{}; // beste treffer zurückgeben
    
    // auto results = seqan3::search(read, jst_adapter{jst_data}, config);
    spm::horspool_matcher matcher{query};

    auto search_tree = libjst::make_volatile(jst_data) | libjst::labelled()
    | libjst::coloured()
    | libjst::trim(spm::window_size(matcher) - 1)
    | libjst::prune()
    | libjst::left_extend(spm::window_size(matcher) - 1)
    | libjst::merge() // make big nodes
    | libjst::seek();

    libjst::tree_traverser_base oblivious_path{search_tree};
    for (auto it = oblivious_path.begin(); it != oblivious_path.end(); ++it) {
         auto && cargo = *it;
         matcher(cargo.sequence(), [&] ([[maybe_unused]] auto && label_finder) {
            std::cout<<"Found hit. Yayy." << cargo.position() <<"\n";
            // callback(query, match_position{.tree_position{cargo.position()},
                                        //    .label_offset{std::ranges::ssize(cargo.sequence()) - seqan2::endPosition(label_finder)}});
        });
}

    std::vector<size_t> hits;
    for (auto&& result : results)
        hits.push_back(result.reference_begin_position());
    
    return hits;
}

void map_reads(std::filesystem::path const & query_path,
    std::filesystem::path const & sam_path,
    const JST_Data & jst_data,
    uint8_t const errors)
{
    sequence_file_t query_file_in{query_path};
    seqan3::sam_file_output sam_out{sam_path, seqan3::fields<seqan3::field::seq,
                                              seqan3::field::id,
                                              seqan3::field::ref_id, 
                                              seqan3::field::ref_offset, 
                                              seqan3::field::cigar, 
                                              seqan3::field::qual, 
                                              seqan3::field::mapq>{}};



    // definiert die Ausrichtungskonfiguration mit einem globalen Alignment-Algorithmus, keine gap penalties
    seqan3::configuration const align_config =
        seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                         seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}} | 
        seqan3::align_cfg::edit_scheme | 
        seqan3::align_cfg::output_alignment{} |
        seqan3::align_cfg::output_begin_position{} | 
        seqan3::align_cfg::output_score{};
 
    for (auto && record : query_file_in) // beginnt Schleife, die jeden Eintrag in der Query-Datei durchläuft
    {
        auto & query = record.sequence(); // referenziert die Sequenz des aktuellen Eintrags
        auto hits = jst_search(jst_data, query);
       
        for (auto hit_pos : hits)
            {
                std::span text_view{std::data(jst_data.sequences[hit_pos]), query.size() +1};

                for (auto&& alignment : seqan3::align_pairwise(std::tie(text_view,query), align_config))
            {
            auto cigar = seqan3::cigar_from_alignment(alignment.alignment()); // wandelt das Alignment in einen CIGAR-String um
            size_t ref_offset = alignment.sequence1_begin_position() + 2; // berechnet Offset in der Referenzsequenz basierend auf der Ausrichtung
            size_t map_qual = 60u + alignment.score(); // berechnet die Mapping-Qualität basierend auf dem Ausrichtungsscore
                //fängt einen neuen Eintrag zur SAM-Datei hinzu
                sam_out.emplace_back(query,
                                     record.id(),
                                     jst_data.ids[hit_pos], 
                                     ref_offset,
                                     cigar,
                                     record.base_qualities(),
                                     map_qual);
            }
        }
    }
}
// Funktion, die das Programm ausfühhrt
void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & query_path,
                 std::filesystem::path const & sam_path,
                 uint8_t const errors)
{
    //JST Ladefunktion???
 auto jst_data = loadjst(reference_path);

map_reads(query_path,
          sam_path, 
          jst_data,
          errors); 
}

// Funktion, die den Argument-Parser konfiguriert
void initialise_argument_parser(seqan3::argument_parser & parser, configuration & args)
{
    parser.info.author = "Julia Bodnar"; // setzt den Autor des Programms
    parser.info.short_description = "JST-based read mapper"; // setzt Beschreibung des Programms 
    parser.info.version = "1.0.0"; // setzt Version des Programms 

    parser.add_option(args.reference_path, // fügt eine Option für den Pfad zur Referenzdatei hinzu
                      'r',
                      "reference",
                      "referenzgenome (JST-format)",
                      sharg::config{.required = true, .validator = sharg::input_file_validator{{"jst"}}});
                    
    parser.add_option(args.query_path, 
                      'q',
                      "query",
                      "Input reads",
                      sharg::config{.required = true, .validator = sharg::input_file_validator{{"fq","fastq"}}});

    parser.add_option(args.sam_path, // fügt eine Option für den Pfad zur SAM Ausgabedatei hinzu
                      'o',
                      "output",
                      "The output SAM file.",
                      sharg::config{.validator = sharg::output_file_validator{{"sam"}}});

    parser.add_option(args.errors, // fügt eine Option für die maximale Anzahl von Feldern hinzu
                      'e',
                      "error",
                      "Maximum allowed errors.",
                      sharg::config{.validator = sharg::arithmetic_range_validator{0, 4}});                 
}

int main(int argc, char const** argv)
{
    sharg::parser parser{"PanMap", argc, argv};
    configuration args{};
    
    initialise_argument_parser(parser, args);
    
    try {
        parser.parse();
    } catch (const sharg::parser_error& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    run_program(args.reference_path, 
               args.query_path, 
               args.sam_path, 
               args.errors);
    
    return 0;
}
