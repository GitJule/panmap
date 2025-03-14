#include <gtest/gtest.h>
#include <loadjst.hpp>

TEST(test, loadrcsstore){
    std::filesystem::path rcsstorepath{std::string{DATADIR}+"local_5refs.jst"};
    rcs_store_t store = loadjst(rcsstorepath);
    EXPECT_EQ(store.size(), 5);
}