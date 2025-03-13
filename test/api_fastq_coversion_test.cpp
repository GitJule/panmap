// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>


#include "app_test.hpp"

// To prevent issues when running multiple API tests in parallel, give each API test unique names:
struct fastq_to_fasta : public app_test
{};

TEST_F(fastq_to_fasta, out_empty)
{GTEST_SKIP();}