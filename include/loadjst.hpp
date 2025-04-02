// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>
#include <cereal/archives/binary.hpp>
#include <fstream>
#include <globaltypes.hpp> 

rcs_store_t loadrcsstore(std::filesystem::path const & rcsstore_path)
{
    std::fstream rcsstream{rcsstore_path};
    if (!rcsstream.good())
      throw std::runtime_error{"Couldn`t read file."};

    rcs_store_t rcsstore{};
    {
        cereal::BinaryInputArchive rcsarchive{rcsstream};
        rcsstore.load(rcsarchive);
    }
    return rcsstore;
}