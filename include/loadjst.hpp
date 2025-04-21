// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cerrno>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>

#include <globaltypes.hpp>

struct JST_Data {
    std::vector<std::string> ids;
    std::vector<seqan3::dna5_vector> sequences;
};

JST_Data loadjst(std::filesystem::path const & path)
{
    JST_Data data;
    seqan3::sequence_file_input fin{path};

    for (auto&&record : fin){
        data.ids.push_back(std::move(record.id()));
        data.sequences.push_back(record.sequence());
    }

    return data;
}

#include <cereal/archives/binary.hpp>

rcs_store_t loadrcsstore(std::filesystem::path const & rcsstore_path)
{
    using namespace std::literals;
    std::ifstream rcsstream{rcsstore_path};
    if (!rcsstream.good())
    {
        throw std::runtime_error{
            "Couldn`t read file: "s + rcsstore_path.string() +
            " (Error: "s + std::strerror(errno) + ")"s};
    }

    rcs_store_t rcsstore{};
    {
        cereal::BinaryInputArchive rcsarchive{rcsstream};
        rcsstore.load(rcsarchive);
    }
    return rcsstore;
}