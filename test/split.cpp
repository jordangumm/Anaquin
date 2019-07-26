#ifdef UNIT_TEST

#include <catch2/catch.hpp>
#include "data/split.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

TEST_CASE("writeSTable_1")
{
    SOptions o;
    o.writer = std::shared_ptr<Writer<>>(new FileWriter("/tmp"));
    writeSTable("test/test3.tsv", "B.tsv", o, 6, 6, 6, "MIX", "TPM");
}

#endif
