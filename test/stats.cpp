#ifdef UNIT_TEST

#include <fstream>
#include <streambuf>
#include <catch2/catch.hpp>
#include "stats/stats.hpp"
#include "tools/tools.hpp"
#include "stats/ss/stats.hpp"

using namespace Anaquin;

TEST_CASE("RAggregateMean")
{
    RAggregateMean("test/test2.tsv", "/tmp/A.tsv", "Mix");
    REQUIRE(readFile("/tmp/A.tsv") == "Mix\tRead\n1\t0\n2\t0.4\n4\t0.4\n8\t4.4\n16\t5.8\n32\t16.4\n64\t27.2\n128\t75.8333\n256\t138\n512\t276.667\n1024\t590.667\n2048\t1312\n4096\t1857.6\n8192\t4044.17\n16384\t6696.2\n32768\t31081.2");
}

TEST_CASE("RGene")
{
    RFilterC("test/rna_sequin.tsv", "/tmp/A1.tsv", std::set<Label> { "GENE", "MIX", "TPM" }, true);
    RAggregateSum("/tmp/A1.tsv", "/tmp/A2.tsv", "GENE", Imputation::ToZero);

    const auto x = readFile("/tmp/A2.tsv");
    
    REQUIRE(!isSubstr(x, "R2_71_"));
    REQUIRE(!isSubstr(x, "R2_72_"));
    REQUIRE(!isSubstr(x, "R2_73_"));
    REQUIRE(!isSubstr(x, "R2_76_"));

    REQUIRE(isSubstr(x, "R2_71"));
    REQUIRE(isSubstr(x, "R2_72"));
    REQUIRE(isSubstr(x, "R2_73"));
    REQUIRE(isSubstr(x, "R2_76"));

    REQUIRE(isSubstr(x, "R1_101\t15.1062\t12033.8"));
    REQUIRE(isSubstr(x, "R1_102\t15.1062\t12311"));
    REQUIRE(isSubstr(x, "R1_103\t966.797\t12957.2"));
}

TEST_CASE("Median")
{
    REQUIRE(med(std::vector<double> { 5079, 4561, 3372, 3371, 4827, 3055 }) == 3966.5);
}

TEST_CASE("RLinear_2")
{
    RGrep("test/meta_sequin_2.tsv", "/tmp/A.txt", "NAME", "MQ_");
    const auto l2 = RLinear("/tmp/A.txt", "NAME", "MIX", "TPM").linear();
    REQUIRE(!std::isnan(l2.p));
    REQUIRE(l2.r == Approx(0.853032));
}

TEST_CASE("Q75")
{
    // Different implementation to R
    REQUIRE(quant(std::vector<double> { 73, 22, 171, 18, 73, 98 }, 0.75) == 85.5);
}

TEST_CASE("RBinaryTSV")
{
    auto r2 = RBinaryTSV("test/A2.tsv", "Mix", "Read");
    auto r3 = RBinaryTSV("test/A3.tsv", "Mix", "Read");

    REQUIRE(r2.size() == 16);
    REQUIRE(r3.size() == 16);

    REQUIRE(r2[1] == 0);
    REQUIRE(r2[2] == 0.4);
    REQUIRE(r2[4] == 0.4);
    REQUIRE(r2[8] == 4.4);
    REQUIRE(r2[16] == 5.8);
    REQUIRE(r2[32] == 16.4);
    REQUIRE(r2[64] == 27.2);
    REQUIRE(r2[128] == 75.8333);
    REQUIRE(r2[256] == 138);

    REQUIRE(r3[1] == 0);
    REQUIRE(r3[2] == 0.894427);
    REQUIRE(r3[4] == 0.894427);
    REQUIRE(r3[8] == 2.60768);
    REQUIRE(r3[16] == 6.18061);
    REQUIRE(r3[32] == 19.4628);
    REQUIRE(r3[64] == 12.1943);
    REQUIRE(r3[128] == 56.226);
    REQUIRE(r3[256] == 66.7742);
}

TEST_CASE("RLinear")
{
    auto s = RLinear("test/mtcars.tsv", "Name", "cyl", "am");
    auto l = s.linear(false, false);
    
    REQUIRE(l.r == Approx(-0.522607));
    REQUIRE(l.R2 == Approx(0.2731181255));
    REQUIRE(l.p == Approx(0.0021512069));
}

TEST_CASE("RFilterR")
{
    RFilterR("test/test.tsv", "/tmp/A.tsv", std::set<std::string> { "Column_1" });
    REQUIRE(readFile("/tmp/A.tsv") == "Column_1\tColumn_2\tColumn_3\nA1\tA2\tA3\nB1\t-\tB3\nC1\tC2\tC3");
    
    RFilterR("test/test.tsv", "/tmp/A.tsv", std::set<std::string> { "Column_2" });
    REQUIRE(readFile("/tmp/A.tsv") == "Column_1\tColumn_2\tColumn_3\nA1\tA2\tA3\nC1\tC2\tC3");

    RFilterR("test/test.tsv", "/tmp/A.tsv", std::set<std::string> { "Column_3" });
    REQUIRE(readFile("/tmp/A.tsv") == "Column_1\tColumn_2\tColumn_3\nA1\tA2\tA3\nB1\t-\tB3\nC1\tC2\tC3");
}

TEST_CASE("RFilterC")
{
    RFilterC("test/test.tsv", "/tmp/A.tsv", std::set<std::string> { "Column_1" });
    REQUIRE(readFile("/tmp/A.tsv") == "Column_2\tColumn_3\nA2\tA3\n-\tB3\n-\t-\nC2\tC3");
    
    RFilterC("test/test.tsv", "/tmp/A.tsv", std::set<std::string> { "Column_2" });
    REQUIRE(readFile("/tmp/A.tsv") == "Column_1\tColumn_3\nA1\tA3\nB1\tB3\n-\t-\nC1\tC3");

    RFilterC("test/test.tsv", "/tmp/A.tsv", std::set<std::string> { "Column_3" });
    REQUIRE(readFile("/tmp/A.tsv") == "Column_1\tColumn_2\nA1\tA2\nB1\t-\n-\t-\nC1\tC2");
}

TEST_CASE("RAggregateSD")
{
    RAggregateSD("test/test2.tsv", "/tmp/A.tsv", "Mix");
    REQUIRE(readFile("/tmp/A.tsv") == "Mix\tRead\n1\t0\n2\t0.894427\n4\t0.894427\n8\t2.60768\n16\t6.18061\n32\t19.4628\n64\t12.1943\n128\t56.226\n256\t66.7742\n512\t78.2628\n1024\t280.816\n2048\t711.865\n4096\t1483.12\n8192\t875.698\n16384\t3273.1\n32768\t10966.7");
}

TEST_CASE("RAggregateMin")
{
    RAggregateMin("test/test2.tsv", "/tmp/A.tsv", "Mix");
    REQUIRE(readFile("/tmp/A.tsv") == "Mix\tRead\n1\t0\n2\t0\n4\t0\n8\t2\n16\t0\n32\t2\n64\t12\n128\t18\n256\t58\n512\t175\n1024\t283\n2048\t484\n4096\t703\n8192\t3055\n16384\t3442\n32768\t15618");
}

TEST_CASE("RAggregateQ50")
{
    RAggregateQ50("test/test2.tsv", "/tmp/A.tsv", "Mix");
    REQUIRE(readFile("/tmp/A.tsv") == "Mix\tRead\n1\t0\n2\t0\n4\t0\n8\t4\n16\t5\n32\t12\n64\t35\n128\t73\n256\t148\n512\t305\n1024\t566.5\n2048\t1382.5\n4096\t1097\n8192\t3966.5\n16384\t5725\n32768\t29790");
}

TEST_CASE("RAggregateQ25")
{
    RAggregateQ25("test/test2.tsv", "/tmp/A.tsv", "Mix");
    
    // Different implementation to R
    REQUIRE(readFile("/tmp/A.tsv") == "Mix\tRead\n1\t0\n2\t0\n4\t0\n8\t2\n16\t2\n32\t4\n64\t16\n128\t20\n256\t65\n512\t180.5\n1024\t319\n2048\t513.5\n4096\t854\n8192\t3213\n16384\t5437\n32768\t26804");
}

TEST_CASE("RAggregateQ75")
{
    RAggregateQ75("test/test2.tsv", "/tmp/A.tsv", "Mix");
    
    // Different implementation to R
    REQUIRE(readFile("/tmp/A.tsv") == "Mix\tRead\n1\t0\n2\t0\n4\t0\n8\t6\n16\t6\n32\t14\n64\t35\n128\t85.5\n256\t189.5\n512\t325\n1024\t718\n2048\t1610.5\n4096\t2417\n8192\t4694\n16384\t6728\n32768\t40950");
}

TEST_CASE("RAggregateMax")
{
    RAggregateMax("test/test2.tsv", "/tmp/A.tsv", "Mix");
    REQUIRE(readFile("/tmp/A.tsv") == "Mix\tRead\n1\t0\n2\t2\n4\t2\n8\t8\n16\t16\n32\t50\n64\t38\n128\t171\n256\t211\n512\t362\n1024\t1025\n2048\t2332\n4096\t4217\n8192\t5079\n16384\t12149\n32768\t42244");
}

TEST_CASE("RMeanCV")
{
    RMeanCV("test/meta_sequin.tsv", "/tmp/A.tsv", "MIX", "READ");
    REQUIRE(readFile("/tmp/A.tsv") == "NAME\tCOUNT\tMEAN\tCV\n1\t0\t0.00\tNA\n2\t1\t0.40\t2.2361\n4\t1\t0.40\t2.2361\n8\t5\t4.4\t0.5927\n16\t4\t5.8\t1.0656\n32\t5\t16.4\t1.1868\n64\t5\t27.2\t0.4483\n128\t6\t75.83\t0.7414\n256\t6\t138\t0.4839\n512\t6\t276.67\t0.2829\n1024\t6\t590.67\t0.4754\n2048\t6\t1312\t0.5426\n4096\t5\t1857.6\t0.7984\n8192\t6\t4044.17\t0.2165\n16384\t5\t6696.2\t0.4888\n32768\t5\t31081.2\t0.3528");
}

TEST_CASE("RLadTable_1")
{
    RLadTable("test/meta_sequin.tsv", "/tmp/A.tsv", "NAME");
    REQUIRE(readFile("/tmp/A.tsv") == "COPY\tMEAN\tSD\tCV\tQ0\tQ25\tQ50\tQ75\tQ100\tRATIO\n1\t87.4\t23.684\t27.0984\t48\t67\t85.5\t107\t119\tNA\n2\t189.5\t53.1439\t28.0443\t117\t140\t176\t230\t263\t2.0585\n4\t347.4\t80.4987\t23.1718\t203\t300.5\t334\t379.5\t502\t1.8977\n8\t668\t166.861\t24.9792\t394\t546\t638.5\t808\t919\t1.9117");
}

#endif
