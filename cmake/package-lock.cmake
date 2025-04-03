# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock (https://github.com/cpm-cmake/CPM.cmake)
# This file should be committed to version control

# cmake-format: off

# sharg
set (PANMAP_SHARG_VERSION c81c1f858054c7114d4d0e82c1c5c2d78574cb5e CACHE STRING "" FORCE)
CPMDeclarePackage (sharg
                   NAME sharg
                   GIT_TAG ${PANMAP_SHARG_VERSION} # main
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "SHARG_NO_TDL ON"
)

# seqan3
set (PANMAP_SEQAN3_VERSION 2863cbbe336a51c21932c69635e814b6e3a8a4ce CACHE STRING "" FORCE)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${PANMAP_SEQAN3_VERSION} # main
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

set (PANMAP_LIBJST_VERSION bb3da768840026bf62a5d9869c04da6c8460a34f)
CPMDeclarePackage (libjst
                   NAME libjst
                   GIT_TAG ${PANMAP_LIBJST_VERSION}
                   GITHUB_REPOSITORY rrahn/libjst
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "LIBJST_DEVELOPER_MODE OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# libspm
set (PAN_MAP_LIBSPM_VERSION 9b1468df9b10913b9e16dd07062d85ba6f948634)
CPMDeclarePackage (libspm
                   NAME libspm
                   GIT_TAG ${PAN_MAP_LIBSPM_VERSION}
                   GITHUB_REPOSITORY rrahn/libspm
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "LIBSPM_DEVELOPER_MODE OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# googletest
set (PANMAP_GOOGLETEST_VERSION 1.16.0 CACHE STRING "" FORCE)
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${PANMAP_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# use_ccache
set (PANMAP_USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37 CACHE STRING "" FORCE)
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${PANMAP_USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)

# cmake-format: on
