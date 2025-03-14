// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <filesystem>

struct configuration
{
    std::filesystem::path reference_path{}; // definiert ein Feld für den Pfad zur Referenzdatei
    std::filesystem::path query_path{}; // für query datei 
    std::filesystem::path index_path{}; // für Index Datei
    std::filesystem::path sam_path{"out.sam"}; // definiert Feld für den Pfad zur SAM-Ausgabedatei mit einem Standardwert
    uint8_t errors{0}; 
    bool verbose{}; // Default is false.
};
