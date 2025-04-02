# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# This includes `cmake/test/declare_datasource.cmake`, which provides the `declare_datasource` function.
include (test/declare_datasource)

# Makes all files in this directory and subdirectories available to the build, i.e., copying them to <build>/data/.
# `datasources.cmake`, `README.md`, and files ending in `.license` are ignored.
# You may organise your data in subdirectories, but each file must have a unique name.
include (test/add_local_data)

# Downloads the file to <build>/data/downloaded.fasta
# The checksum of the downloaded file is checked to ensure that:
# * The download was not corrupted.
# * The file was not tampered with.
#
# Some ways to obtain the checksum:
#
#   wget --quiet -O- https://ftp.seqan.de/app-template/downloaded.fasta | sha256sum
#   c3cb990ca1a25c7e31be3c6c2d009238d9ac9a44b2b7c143753c1e2881699077  -
#   ^------------------------- checksum ---------------------------^  ^ file (- is stdin, the file was piped)
#
# You can also use curl:
#   curl --silent https://ftp.seqan.de/app-template/downloaded.fasta | sha256sum
#
# If you have the file locally, you can use `sha256sum` directly:
#   sha256sum downloaded.fasta
#   c3cb990ca1a25c7e31be3c6c2d009238d9ac9a44b2b7c143753c1e2881699077  downloaded.fasta
declare_datasource (FILE downloaded.fasta # This is a custom name. It does not have to match the file name.
                    URL https://ftp.seqan.de/app-template/downloaded.fasta
                    URL_HASH SHA256=c3cb990ca1a25c7e31be3c6c2d009238d9ac9a44b2b7c143753c1e2881699077
)

declare_datasource (FILE local_5refs.jst # This is a custom name. It does not have to match the file name.
                    URL ${CMAKE_SOURCE_DIR}/test/data/test.jst
                   # URL_HASH SHA256=ce5b3c884df5eaa39326b1270250946185422071b2f6d79f137ea184a2f66c32
             
                   URL_HASH SHA256=b0f8aff1776fcd9106eee7c135573ee0322e877768b94afb1815d539fb454a3c
)
