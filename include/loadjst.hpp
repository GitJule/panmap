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

