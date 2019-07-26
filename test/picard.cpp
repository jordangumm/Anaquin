#include <catch2/catch.hpp>
#include "tools/picard.hpp"
#include "data/standard.hpp"
#include "Genomics/Genomics.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

TEST_CASE("Picard_3")
{
    UserReference r;
    
    r.v4 = std::shared_ptr<VCFLadder>(
            new VCFLadder(Standard::addVCF(Reader("data/genome/chrQ/sequin_smallvariants_chrQS_2.6.vcf"),
                                                                    nullptr)));
    RegionOptions o;
    auto r1 = BedData(Standard::readBED(Reader("data/genome/hg38/sequin_regions_hg38_2.6.bed"), o)).inters();
    auto r2 = BedData(Standard::readBED(Reader("data/genome/chrQ/sequin_regions_chrQS_2.6.bed"), o)).inters();

    Standard::instance().gen.finalize(r);
    assert(Standard::instance().gen.v4());
    
    std::map<ChrID, Sequence> chrs;
    ParserFA::parse(Reader("data/genome/sequin_sequences_2.6.fa"), [&](const ParserFA::Data &x) {
        chrs[x.id] = x.seq;
    });
    
    Picard p(chrs);
    
    PicardOption o2;
    o2.ignoreSkipClip = true;
    std::vector<DInter *> v;
    
    ParserBAM::parse("test/calibrate_calibrated.bam", [&](ParserBAM::Data &x, const ParserBAM::Info &) {
        v.clear();
        if (r1.count(x.cID) && (r1.at(x.cID).contains(x.l, &v)))
        {
            DInter *d = nullptr;
            for (auto &i : v)
            {
                if (!isSubstr(i->id(), "_A"))
                {
                    d = i;
                }
            }
            
            if (!d)
            {
                return;
            }
            
            const auto sID = d->id();

            if (!chrs.count(sID))
            {
                return;
            }
            
            x.lName(); x.lCigar(); x.lSeq();
            
            assert(!d->id().empty());
            assert(chrs.count(sID));
            
            /*
             * Convert genome coordinate into sequin coordinate
             */
            
            auto t = x;
            
            if (x.isReverseC)
            {
                t.l.start = chrs.at(sID).size() - (x.l.start - d->l().start) - t.seq.size() + 1;
                std::reverse(t.seq.begin(), t.seq.end());
            }
            else
            {
                t.l.start = chrs.at(sID).size() - (x.l.start - d->l().start) - t.seq.size() + 1;
                std::reverse(t.seq.begin(), t.seq.end());
            }
            
            if (t.l.start >= 0)
            {
                p.analyze(sID, t, o2);
            }
        }
    });
    
    REQUIRE(p.skips == 0);
    REQUIRE(p.sumIL(1) == 0);
    REQUIRE(p.sumDL(1) == 0);
    REQUIRE(p.sumSB(SNPBin::GA) == 14);
    REQUIRE(p.sumSB(SNPBin::AC) == 12);
    REQUIRE(p.sumSB(SNPBin::GT) == 11);
    REQUIRE(p.sumSB(SNPBin::CT) == 8);
    REQUIRE(p.sumSB(SNPBin::TC) == 14);
    REQUIRE(p.sumSB(SNPBin::Match) == 159);    
}

TEST_CASE("Picard_1")
{
    UserReference r;

    r.v4 = std::shared_ptr<VCFLadder>(
        new VCFLadder(Standard::addVCF(Reader("data/genome/chrQ/sequin_smallvariants_chrQS_2.6.vcf"),
            nullptr)));
 
    Standard::instance().gen.finalize(r);
    assert(Standard::instance().gen.v4());
    
    std::map<ChrID, Sequence> chrs;
    ParserFA::parse(Reader("data/genome/chrQ/genome_chrQ_decoys_2.6.fa"), [&](const ParserFA::Data &x) {
        if (x.id == GENOMICS_DECOY_CHROM) { chrs[x.id] = x.seq; }
    });
    
    Picard p(chrs);
    ParserBAM::parse("test/perfect.bam", [&](ParserBAM::Data &x, const ParserBAM::Info &) {
        x.lName(); x.lCigar(); x.lSeq();
        if (x.cID == GENOMICS_DECOY_CHROM) { p.analyze(GENOMICS_DECOY_CHROM, x); }
    });
    
    REQUIRE(p.skips == 7);
    //REQUIRE(p.sumIL( ins.empty());
    //REQUIRE(p.dels.empty());
    REQUIRE(p.sumSB(SNPBin::GA) == 3);
    REQUIRE(p.sumSB(SNPBin::CT) == 1);
    REQUIRE(p.sumSB(SNPBin::TC) == 1);
    REQUIRE(p.sumSB(SNPBin::AC) == 0);
    REQUIRE(p.sumSB(SNPBin::GT) == 0);
    REQUIRE(p.sumSB(SNPBin::Match) == 121047);    
}

TEST_CASE("Picard_2")
{
    UserReference r;
    
    r.v4 = std::shared_ptr<VCFLadder>(
            new VCFLadder(Standard::addVCF(Reader("data/genome/chrQ/sequin_smallvariants_chrQS_2.6.vcf"),
                nullptr)));
    
    Standard::instance().gen.finalize(r);
    assert(Standard::instance().gen.v4());

    std::map<ChrID, Sequence> chrs;
    ParserFA::parse(Reader("data/genome/chrQ/genome_chrQ_decoys_2.6.fa"), [&](const ParserFA::Data &x) {
        chrs[x.id] = x.seq;
    });

    Picard p(chrs);    
    ParserBAM::parse("test/chrQ.sam", [&](ParserBAM::Data &x, const ParserBAM::Info &) {
        if (x.cID == GENOMICS_DECOY_CHROM)
        {
            x.lName(); x.lCigar(); x.lSeq();
            p.analyze(GENOMICS_DECOY_CHROM, x);
        }
    });
    
    REQUIRE(p.skips == 124);
    REQUIRE(p.sumIL(1) == 1);
    REQUIRE(p.sumDL(1) == 2);
    REQUIRE(p.sumSB(SNPBin::GA) == 1);
    REQUIRE(p.sumSB(SNPBin::AC) == 4);
    REQUIRE(p.sumSB(SNPBin::GT) == 3);
    REQUIRE(p.sumSB(SNPBin::CT) == 2);
    REQUIRE(p.sumSB(SNPBin::TC) == 3);
    REQUIRE(p.sumSB(SNPBin::Match) == 1461);
}
