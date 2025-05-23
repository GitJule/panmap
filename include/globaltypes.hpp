// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <libspm/seqan/alphabet.hpp>
#include <libjst/rcms/compressed_multisequence.hpp> 
#include <libjst/rcms/rcs_store.hpp>
#include <libjst/coverage/int_coverage.hpp> 

using alphabet_t = spm::dna5;
using coverage_t = libjst::int_coverage<uint32_t>;
using reference_t = std::vector<alphabet_t>;

using cms_t = libjst::compressed_multisequence<reference_t, coverage_t>;
using rcs_store_t = libjst::rcs_store<reference_t, cms_t>;
