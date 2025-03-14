// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sharg/all.hpp>

#include <configuration.hpp>
/*
int main(int argc, char ** argv)
{
    // Configuration
    configuration config{};

    // Parser
    sharg::parser parser{"PanMap", argc, argv};

    // General information.
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    // Positional option: The FASTQ file to convert.
    parser.add_positional_option(config.fastq_input,
                                 sharg::config{.description = "The FASTQ file containing the reads",
                                               .validator = sharg::input_file_validator{{"fq", "fastq"}}});

    // Open: Output FASTA file. Default: print to terminal 
    parser.add_option(config.fasta_output,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output",
                                    .description = "The output SAM file.",
                                    .default_message = "Print to terminal (stdout)",
                                    .validator = sharg::output_file_validator{}});

    // Flag: Verose output.
    parser.add_flag(
        config.verbose,
        sharg::config{.short_id = 'v', .long_id = "verbose", .description = "Give more detailed information."});

    try
    {
        parser.parse(); // Trigger command line parsing.
    }
    catch (sharg::parser_error const & ext) // Catch user errors.
    {
        std::cerr << "Parsing error. " << ext.what() << '\n'; // Give error message.
        return -1;
    }


    if (config.verbose) // If flag is set.
        std::cerr << "Mapping was successful. Congrats!\n";

    return 0;
}

*/

#    include <fstream>
#    include <span>
 
#    include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#    include <seqan3/alignment/configuration/all.hpp>
#    include <seqan3/alignment/pairwise/align_pairwise.hpp>
#    include <seqan3/argument_parser/all.hpp>
#    include <seqan3/io/sam_file/output.hpp>
#    include <seqan3/io/sequence_file/input.hpp>
#    include <seqan3/search/all.hpp>
#    include <seqan3/search/fm_index/bi_fm_index.hpp>
 
#    include <cereal/archives/binary.hpp>
 
struct reference_storage_t
{   std::vector<std::string> ids;// Vektor zur Speicherung der Identifikatoren der Referenzsequenzen
    std::vector<std::vector<seqan3::dna5>> seqs; // Vektor von Vektoren zur Speicherung der Referenzsequenzen selbst
};
 
void read_reference(std::filesystem::path const & reference_path, reference_storage_t & storage) // Funktion, die Referenzdatei liest und Daten in reference_storage_t-Struktur speichert
{
    seqan3::sequence_file_input reference_in{reference_path}; // öffnet Referenzdatei
    for (auto && record : reference_in) // beginnt eine Schleife, die jeden Eintrag in der Referenzdatei durchläuft
    {
        storage.ids.push_back(record.id()); // fügt die Identifikator des aktuellen Eintrags zur Liste der Identifikatoren hinzu
        storage.seqs.push_back(record.sequence()); // fügt die Sequenz des aktuellen Eintrags zur Liste der Sequenzen hinzu
    }
}

// Funktion, die Reads gegen eine Referenz mappt
 
void map_reads(std::filesystem::path const & query_path, // Funktion, die Reads gegen eine Referenz mappt
               std::filesystem::path const & index_path, 
               std::filesystem::path const & sam_path,
               reference_storage_t & storage,
               uint8_t const errors)
{
    seqan3::bi_fm_index<seqan3::dna5, seqan3::text_layout::collection> index; // deklariert einen bidirektionalen FM-Index für DNA-Sequenzen
    {
        std::ifstream is{index_path, std::ios::binary}; //  öffnet die Indexdatei im binären Modus
        cereal::BinaryInputArchive iarchive{is}; // erstellt einen binären Archiv-Reader für die Serialisierung
        iarchive(index); // lädt den Index aus der Datei
    }
 
    seqan3::sequence_file_input query_file_in{query_path}; //  öffnet die Query-Datei
 
    //  initialisiert die SAM-Ausgabedatei mit bestimmten Feldern
    seqan3::sam_file_output sam_out{sam_path,
                                    seqan3::fields<seqan3::field::seq,
                                                   seqan3::field::id,
                                                   seqan3::field::ref_id,
                                                   seqan3::field::ref_offset,
                                                   seqan3::field::cigar,
                                                   seqan3::field::qual,
                                                   seqan3::field::mapq>{}};

    // definiert die Suchkonfiguration mit einer maximalen Anzahl von Fehlern und der Option, alle besten Treffer zu berücksichtigen
    seqan3::configuration const search_config =
        seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{errors}}
        | seqan3::search_cfg::hit_all_best{};

    // definiert die Ausrichtungskonfiguration mit einem globalen Alignment-Algorithmus und bestimmten Optionen zur Ausgabe
    seqan3::configuration const align_config =
        seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                         seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
        | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_alignment{}
        | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_score{};
 
    for (auto && record : query_file_in) // beginnt Schleife, die jeden Eintrag in der Query-Datei durchläuft
    {
        auto & query = record.sequence(); // referenziert die Sequenz des aktuellen Eintrags
        for (auto && result : search(query, index, search_config)) // beginnt Schleife, die die Suchergebnisse für die aktuelle Sequenz durchläuft
        {
            size_t start = result.reference_begin_position() ? result.reference_begin_position() - 1 : 0; // berechnet Startpunkt der Referenzsequenz basierend auf dem Suchergebnis
            std::span text_view{std::data(storage.seqs[result.reference_id()]) + start, query.size() + 1}; // erstellt Span, der Teil der Referenzsequenz darstellt, beginnend am berechneten Startpunkt
 
            for (auto && alignment : seqan3::align_pairwise(std::tie(text_view, query), align_config)) // beginnt Schleife, die Ausrichtungsergebnisse für die aktuelle Sequenz durchläuft
            {
                auto cigar = seqan3::cigar_from_alignment(alignment.alignment()); // wandelt das Alignment in einen CIGAR-String um
                size_t ref_offset = alignment.sequence1_begin_position() + 2 + start; // berechnet Offset in der Referenzsequenz basierend auf der Ausrichtung
                size_t map_qual = 60u + alignment.score(); // berechnet die Mapping-Qualität basierend auf dem Ausrichtungsscore
                //fügt einen neuen Eintrag zur SAM-Datei hinzu
                sam_out.emplace_back(query,
                                     record.id(),
                                     storage.ids[result.reference_id()],
                                     ref_offset,
                                     cigar,
                                     record.base_qualities(),
                                     map_qual);
            }
        }
    }
}
// Funktion, die das Programm ausführt
void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & query_path,
                 std::filesystem::path const & index_path,
                 std::filesystem::path const & sam_path,
                 uint8_t const errors)
{
    reference_storage_t storage{}; //  initialisiert die Referenzspeicherung
    read_reference(reference_path, storage); // liest die Referenzdatei.
    map_reads(query_path, index_path, sam_path, storage, errors); // mappt die Reads gegen die Referenz
}
 
struct cmd_arguments
{
    std::filesystem::path reference_path{}; // definiert ein Feld für den Pfad zur Referenzdatei
    std::filesystem::path query_path{}; // für query datei 
    std::filesystem::path index_path{}; // für Index Datei
    std::filesystem::path sam_path{"out.sam"}; // definiert Feld für den Pfad zur SAM-Ausgabedatei mit einem Standardwert
    uint8_t errors{0}; // definiert ein Feld für die maximale Anzahl von Fehlern
};
// Funktion zur Initialisierung des Argument-Parsers
// Funktion, die den Argument-Parser konfiguriert
void initialise_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "E. coli"; // setzt den Autor des Programms
    parser.info.short_description = "Map reads against a reference."; // setzt Beschreibung des Programms 
    parser.info.version = "1.0.0"; // setzt Version des Programms 
    parser.add_option(args.reference_path, // fügt eine Option für den Pfad zur Referenzdatei hinzu
                      'r',
                      "reference",
                      "The path to the reference.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta"}});
    parser.add_option(args.query_path, // fügt eine Option für den Pfad zur Query Datei hinzu
                      'q',
                      "query",
                      "The path to the query.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fq", "fastq"}});
    parser.add_option(args.index_path, // fügt eine Option für den Pfad zur Indexdatei hinzu
                      'i',
                      "index",
                      "The path to the index.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"index"}});
    parser.add_option(args.sam_path, // fügt eine Option für den Pfad zur SAM Ausgabedatei hinzu
                      'o',
                      "output",
                      "The output SAM file path.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"sam"}});
    parser.add_option(args.errors, // fügt eine Option für die maximale Anzahl von Feldern hinzu
                      'e',
                      "error",
                      "Maximum allowed errors.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 4});
}
// main function
int main(int argc, char const ** argv)
{
    seqan3::argument_parser parser("Mapper", argc, argv);
    cmd_arguments args{};
 
    initialise_argument_parser(parser, args);
 
    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
        return -1;
    }
 
    run_program(args.reference_path, args.query_path, args.index_path, args.sam_path, args.errors);
 
    return 0;
}
