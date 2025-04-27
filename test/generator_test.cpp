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

TEST(test, generator) {
    auto store = loadjst("test.jst");
    
    generate_test_sequences( 
        store,
        "test_queries.fasta",
        50,   // Patternlänge
        1000, // Samples pro Knoten
        1234  // Fester Seed
    );
    
    seqan3::sequence_file_input fasta_in{"test_queries.fasta"};
    size_t count = 0;
    for (auto const& record : fasta_in) {
        EXPECT_EQ(record.sequence().size(), 50); 
        count++;
    }
    EXPECT_GE(count, 1000);
}