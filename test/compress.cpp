#include <fstream>
#include <catch.hpp>
#include "data/reader.hpp"
#include "tools/tools.hpp"

using namespace Anaquin;

TEST_CASE("Compress_1")
{
    Reader r("test/hello.txt.gz");
    std::string s1, s2, s3;
    
    const auto r1 = r.nextLine(s1);
    const auto r2 = r.nextLine(s2);
    const auto r3 = r.nextLine(s3);

    REQUIRE(r1);
    REQUIRE(r2);
    REQUIRE(!r3);
    REQUIRE(s1 == "Hello World!");
    REQUIRE(s2 == "Hi Again!");
}

TEST_CASE("Compress_2")
{
    system("rm -f /tmp/A1.txt.gz");
    system("rm -f /tmp/A2.txt.gz");
    
    auto w0 = std::shared_ptr<std::ofstream>(new std::ofstream("/tmp/A1.txt.gz", std::ios::binary | std::ios::out));
    auto w1 = std::shared_ptr<std::ofstream>(new std::ofstream("/tmp/A2.txt.gz", std::ios::binary | std::ios::out));
    
    Compressor c0(w0, 5);
    Compressor c1(w1, 5);
    
    for (auto i = 0; i < 5; i++) { c0.write("0"); }
    for (auto i = 0; i < 5; i++) { c1.write("1"); }
    
    w0->close(); w1->close();
    
    system("gunzip -f /tmp/A1.txt.gz");
    system("gunzip -f /tmp/A2.txt.gz");
    
    std::ifstream f1("/tmp/A1.txt");
    std::ifstream f2("/tmp/A2.txt");
    
    std::string x1((std::istreambuf_iterator<char>(f1)), std::istreambuf_iterator<char>());
    std::string x2((std::istreambuf_iterator<char>(f2)), std::istreambuf_iterator<char>());
    
    REQUIRE(x1 == "00000");
    REQUIRE(x2 == "11111");
}

TEST_CASE("Compress_3")
{
    const auto x1 = "ABCD1234";
    
    system("rm -f /tmp/A.txt");
    std::ofstream f1;
    f1.open("/tmp/A.txt.gz");
    f1 << compressGZ(x1);
    f1.close();
    system("gunzip -f /tmp/A.txt.gz");
    
    std::ifstream f2("/tmp/A.txt");
    std::string x2((std::istreambuf_iterator<char>(f2)), std::istreambuf_iterator<char>());
    REQUIRE(x1 == x2);
}
