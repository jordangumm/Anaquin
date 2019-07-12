#include <catch.hpp>
#include "tools/tools.hpp"

using namespace Anaquin;

TEST_CASE("Tool_1")
{
    //REQUIRE(!headFQGZ("test/partition_B8_10_1.fq.gz.tmp"));
    //REQUIRE(!headFQGZ("test/partition_B8_10_2.fq.gz.tmp"));
    REQUIRE(headFQGZ("test/mini.1.fastq.gz"));
    REQUIRE(headFQGZ("test/mini.2.fastq.gz"));
    REQUIRE(headFQGZ("test/partition_B12_12_1.fq.gz.tmp"));
}

TEST_CASE("Tool_2")
{
    REQUIRE(!headFQGZ("test/partition_B4_10_1.fq.gz.tmp"));
    const auto x = hackUZero("test/partition_B4_10_1.fq.gz.tmp");
    REQUIRE(!headFQGZ("test/partition_B4_10_1.fq.gz.tmp"));
    REQUIRE(headFQGZ(x));
}
