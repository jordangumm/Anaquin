#include "tools/report.hpp"
#include "Genomics/g_split.hpp"
#include "Genomics/Genomics.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

typedef GSplit::Stats Stats;
typedef GSplit::Options Options;

static std::string bin2Label(Bin x)
{
    switch (x)
    {
        case MI: { return "MSI";        }
        case HP: { return "HP";         }
        case IF: { return "Info";       }
        case MT: { return "Mito";       }
        case HL: { return "HLA";        }
        case LD: { return "Ladder";     }
        case IM: { return "Immune";     }
        case ES: { return "Sample";     }
        case VC: { return "Vector";     }
        case SV: { return "Structural"; }
        case GR: { return "Germline";   }
        case SO: { return "Somatic";    }
        case MS: { return "Micro";      }
        default: { throw std::runtime_error("Unknown: " + std::to_string(x)); }
    }
}

/*
 * Building allele frequency ladder. Unlike abundance quantification, this must be
 * done by unique k-mers.
 */

void GSplit::buildAF(Stats &stats, const Options &o)
{
    o.info("Building allele frequency ladder");
    const auto x = Standard::instance().gen.v1()->data.seq2Var();

    for (const auto &std : stats.K.stds)
    {
        if (GBin(std) != SO)
        {
            continue;
        }

        for (const auto &seq : stats.K.aSeqs[std])
        {
            if (x.count(__hack__(seq)))
            {
                const auto R = stats.K.rSeqs[std];
                const auto V = seq;
                
                const auto hasR = stats.K.uniqs.count(R) && stats.K.uniqs.at(R).size() >= 3;
                const auto hasV = stats.K.uniqs.count(V) && stats.K.uniqs.at(V).size() >= 3;
                
                stats.VR[seq] = hasR ? med(toVector(stats.K.uniqs.at(R))) : NAN;
                stats.VV[seq] = hasV ? med(toVector(stats.K.uniqs.at(V))) : NAN;
            }
        }
    }

    const auto &r = Standard::instance().gen;
    const auto seq2Var = Standard::instance().gen.v1()->data.seq2Var();
    
    for (const auto &std : stats.K.stds)
    {
        if (GBin(std) != LD)
        {
            const auto R = stats.K.rSeqs.count(std) ? stats.K.rSeqs[std] : "";
            
            // There could be more than a single variant sequin
            for (const auto &seq : stats.K.aSeqs[std])
            {
                const auto V = seq;
                
                // Median k-mer for reference sequin
                const auto rn = stats.K.uniqs.count(R) ? med(toVector(stats.K.uniqs.at(R))) : 0;
                
                // Median k-mer for variant sequin
                const auto vn = stats.K.uniqs.count(V) ? med(toVector(stats.K.uniqs.at(V))) : 0;
                
                if (!r.l1()->contains(__hack__(V), Mix_1))
                {
                    continue;
                }
                else if (rn + vn == 0)
                {
                    if (!R.empty()) { o.logInfo(R + " not detected"); }
                    if (!V.empty()) { o.logInfo(V + " not detected"); }
                    continue;
                }
                else if (vn == 0)
                {
                    o.logInfo(V + " abundance is zero");
                    continue;
                }
                
                // Expected allele frequency
                const auto exp = r.l1()->input(__hack__(V), Mix_1);
                
                // Measured allele frequency by unique k-mers
                const auto obs = (float) vn / (rn + vn);
                
                if (seq2Var.count(__hack__(V)))
                {
                    if (isSoma(std))      { stats.af.add(seq, exp, obs); }
                    else if (isGerm(std)) { stats.gm.add(seq, exp, obs); }
                }
            }
        }
    }    
}

Stats GSplit::analyze(const FileName &f1, const FileName &f2, const Options &o)
{
    o.info("Index: "     + o.index);
    o.info("Kmer: "      + toString(o.k));
    o.info("Skip: "      + toString(o.skipKM));
    o.info("FASTQ: "     + toString(o.forceFQ));
    o.info("Threads: "   + toString(o.thr));
    o.info("Threshold: " + toString(o.rule));

    Stats stats;    

    // Running Kallisto
    SAlign(f1, stats, o, nullptr, true); SKallisto(stats, f1, f2, o);

    // Build somatic allele frequency ladder
    buildAF(stats, o);

    // Build synthetic ladder
    buildSL(stats, o);
    
    return stats;
}

void GSplit::buildSL(Stats &stats, const Options &)
{
    const auto &r = Standard::instance().gen;
    
    for (const auto &std : stats.K.stds)
    {
        if (GBin(std) == LD)
        {
            if (stats.K.aSeqs.count(std))
            {
                for (const auto &seq : stats.K.aSeqs.at(std))
                {
                    if (r.l3()->contains(seq) && stats.R.d2u.meds.count(seq))
                    {
                        stats.ld.add(seq, r.l3()->input(seq), stats.R.d2u.meds.at(seq));
                    }
                }
            }
        }
    }
}

void GSplit::writeGReport(const FileName &, const SOptions &o)
{
    if (o.report)
    {
        const auto tmp = o.work + "/" + o.name +  "_features.bed";
        copy(Standard::instance().gen.a1()->src, tmp);        
        typedef Report::Options::KMOptions KOptions;
        Report::Options o2(o);
        o2.k = std::shared_ptr<KOptions>(new KOptions());
        Report::genome(o2);
        removeF(tmp);
    }
}

static void writeSummary(const FileName &file, const FileName &f1, const FileName &f2, const Stats &stats, const Options &o)
{
    o.generate(file);
    o.writer->open(file);

    const auto f = "SEQUIN REPORT:                 %1%\n\n"
                   "REFERENCE FILES\n"
                   "Reference index:               %2%\n\n"
                   "LIBRARY INFORMATION\n"
                   "Version:                       %3%\n"
                   "Instrument ID:                 %4%\n"
                   "Run number:                    %5%\n"
                   "Flowcell ID:                   %6%\n"
                   "Lane:                          %7%\n\n"
                   "USER-SUPPLIED FILES\n"
                   "Input file (first):            %8%\n"
                   "Input file (second):           %9%\n\n"
                   "ANAQUIN PARAMETERS\n"
                   "K-mer length:                  %10%\n"
                   "Threshold:                     %11%\n\n"
                   "PARTITION SUMMARY\n"
                   "Sample reads:                  %12% (%13%%%)\n"
                   "Sequin reads:                  %14% (%15%%%)\n"
                   "Ladder reads:                  %16% (%17%%%)\n"
                   "Structural reads:              %18% (%19%%%)\n"
                   "Immune reads:                  %20% (%21%%%)\n"
                   "HLA reads:                     %22% (%23%%%)\n"
                   "Mitochondria reads:            %24% (%25%%%)\n"
                   "Vector reads:                  %26% (%27%%%)\n"
                   "Information reads:             %28% (%29%%%)\n"
                   "Dilution:                      %30%%%\n"
                   "Total reads:                   %31%\n\n"
                   "OUTPUT FILES\n"
                   "Sample reads path:             %32%\n"
                   "Sequin reads path:             %33%\n"
                   "Ladder reads path:             %34%\n"
                   "Structural reads path:         %35%\n"
                   "Immune reads path:             %36%\n"
                   "Mitochondria reads path:       %37%\n"
                   "HLA reads path:                %38%\n"
                   "Information reads path:        %39%\n"
                   "Vector reads path:             %40%\n\n"
                   "LADDER - LIBRARY QUALITY\n"
                   "Slope:                         %41%\n"
                   "R2:                            %42%\n"
                   "Mean Ratio:                    %43%\n"
                   "Ladder table:                  %44%\n\n"
                   "SEQUIN SOMATIC QUANTIFICATION\n"
                   "Slope:                         %45%\n"
                   "R2:                            %46%\n"
                   "Sequin table:                  %47%";

    LinearModel af, ld;
    try { af = stats.af.linear(true, true); } catch(...) {}
    try { ld = stats.ld.linear(true, true); } catch(...) {}
    
    const auto tmp = tmpFile();
    
    // Only ladder sequins
    RGrep(o.work + "/split_ladder.tsv", tmp, "NAME", "LD_"); const auto l1 = RLinear(tmp, "NAME", "UNIT", "READ").linear();
    
    // Only somatic sequins
    RGrep(o.work + "/split_sequin.tsv", tmp, "LABEL", "Somatic"); const auto l2 = RLinear(tmp, "NAME", "EXP_FREQ", "OBS_FREQ").linear();

    const auto &r = Standard::instance().gen;
    
    #define C(x) (stats.K.c1.count(x) ? stats.K.c1.at(x) : 0)
    #define P(x) (100.0 * stats.K.binP(x))

    const auto gn = C(GR) + C(SO) + C(MS) + C(HP);
    const auto gp = P(GR) + P(SO) + P(MS) + P(HP);
    assert(gp >= 0 && gp <= 100);
    
    extern FASTQ __KFQInfo__;
    const auto fo = __KFQInfo__.format();

    o.writer->write((boost::format(f) % date()
                                      % o.index
                                      % SVersion(Standard::instance().gen, stats.K)
                                      % __KFQInfo__.inst(fo) // 4
                                      % __KFQInfo__.run(fo)  // 5
                                      % __KFQInfo__.flow(fo) // 6
                                      % __KFQInfo__.lane(fo) // 7
                                      % f1                   // 8
                                      % f2                   // 9
                                      % o.k                  // 10
                                      % o.rule               // 11
                                      % C(ES)                // 12
                                      % S2(100.0 * stats.K.binP(ES))
                                      % gn                // 14
                                      % S2((Proportion) gn / stats.K.total())
                                      % C(LD)                // 16
                                      % S2(100.0 * stats.K.binP(LD))
                                      % C(SV)                // 18
                                      % S2(100.0 * stats.K.binP(SV))
                                      % C(IM)                // 20
                                      % S2(100.0 * stats.K.binP(IM))
                                      % C(HL)                // 22
                                      % S2(100.0 * stats.K.binP(HL))
                                      % C(MI)                // 24
                                      % S2(100.0 * stats.K.binP(MI))
                                      % C(VC)               // 26
                                      % S2(100.0 * stats.K.binP(VC))
                                      % C(IF)               // 28
                                      % S2(100.0 * stats.K.binP(IF))
                                      % (100.0 * stats.dil())
                                      % stats.K.total()             // 31
                                      % (o.work + "/split_sample*") // 32
                                      % (o.work + "/split_sequin*") // 33
                                      % (o.work + "/split_ladder*") // 34
                                      % (o.work + "/split_sv*")     // 35
                                      % (o.work + "/split_immune*") // 36
                                      % (o.work + "/split_mito*")   // 37
                                      % (o.work + "/split_hla*")    // 38
                                      % (o.work + "/split_info*")   // 39
                                      % (o.work + "/split_vector*") // 40
                                      % replaceNA(l1.m)  // 41
                                      % replaceNA(l1.R2) // 42
                                      % replaceNA(RLadTable(o.work + "/split_ladder.tsv", tmp, "NAME")) // 43
                                      % (o.work + "/split_ladder_table.tsv") // 44
                                      % replaceNA(l2.m)  // 45
                                      % replaceNA(l2.R2) // 46
                                      % (o.work + "/split_sequin_table.tsv") // 47
                     ).str());
    o.writer->close();
}

void GSplit::writeQuin(const FileName &file, const Stats &stats, const SOptions &o)
{
    const auto &r = Standard::instance().gen;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%\t%14%\t%15%\t%16%%17%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "NAME"
                                           % "LABEL"
                                           % "CHROM"
                                           % "POSITION"
                                           % "EXP_FREQ"
                                           % "OBS_FREQ"
                                           % "REF_COUNT"
                                           % "VAR_COUNT"
                                           % "MIN"
                                           % "Q25"
                                           % "Q50"
                                           % "Q75"
                                           % "MAX"
                                           % "READ"
                                           % "GENOTYPE"
                                           % "TYPE"
                                           % r.a1()->strForKeys()).str());
    
    std::map<SequinID, Base> r1;
    ParserFA::size(o.index, r1);

    #define _S1_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? MISSING : toString(x.at(y)) : MISSING)
    #define _S2_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? MISSING : toString(x.at(y), 2) : MISSING)

    for (const auto &std : stats.K.stds)
    {
        if (isVector(std))
        {
            continue;
        }
        
        auto write = [&](const SequinID &x)
        {
            // Can we find it in VCF? Multiple variants possible.
            const auto v = r.v1()->data.find(__hack__(x));

            /*
             * Never give out number of reads for somatic and germline, because a read can span multiple
             * alleles within the same standard. We don't know unique number of reads for each allele.
             */
            
            const auto read = !isSoma(x) && !isGerm(x) ? stats.K.sqc.count(x) ? toString(stats.K.sqc.at(x)) : MISSING : MISSING;

            if (!v)
            {
                o.writer->write((boost::format(format) % x
                                                       % bin2Label(GBin(GSeq2Std(x)))
                                                       % MISSING
                                                       % MISSING
                                                       % MISSING
                                                       % MISSING
                                                       % MISSING
                                                       % MISSING
                                                       % _S1_(stats.R.d2u.mins, x)
                                                       % _S2_(stats.R.d2u.q25,  x)
                                                       % _S2_(stats.R.d2u.meds, x)
                                                       % _S2_(stats.R.d2u.q75,  x)
                                                       % _S1_(stats.R.d2u.maxs, x)
                                                       % read
                                                       % MISSING
                                                       % MISSING
                                                       % r.a2()->strForNulls()).str());
                return;
            }
            
            const auto R = stats.VR.count(x) ? stats.VR.at(x) : NAN; // Counts for reference allele
            const auto V = stats.VV.count(x) ? stats.VV.at(x) : NAN; // Counts for variant allele
            
            // Expected allele frequency
            const auto exp = r.l1()->contains(__hack__(x), Mix_1) ? r.l1()->input(__hack__(x), Mix_1) : NAN;
            
            // Measured allele frequency
            const auto obs = stats.af.count(x) ? stats.af.at(x).y : stats.gm.count(x) ? stats.gm.at(x).y : NAN;
            
            o.writer->write((boost::format(format) % x
                                                   % bin2Label(GBin(GSeq2Std(x)))
                                                   % v->cID
                                                   % v->l.start
                                                   % replaceNA(exp, 6)
                                                   % replaceNA(obs, 6)
                                                   % replaceNA(R)
                                                   % replaceNA(V)
                                                   % _S1_(stats.R.d2u.mins, x)
                                                   % _S2_(stats.R.d2u.q25,  x)
                                                   % _S2_(stats.R.d2u.meds, x)
                                                   % _S2_(stats.R.d2u.q75,  x)
                                                   % _S1_(stats.R.d2u.maxs, x)
                                                   % read
                                                   % gt2str(v->gt)
                                                   % var2str(v->type())
                                                   % r.a2()->strForLocus(v->cID, v->l)).str());
        };

        if (stats.K.aSeqs.count(std))
        {
            for (const auto &seq : stats.K.aSeqs.at(std))
            {
                write(seq);
            }
        }
    }
    
    o.writer->close();
}

void GSplit::report(const FileName &f1, const FileName &f2, const Options &o)
{
    const auto stats = analyze(f1, f2, o);
    
    // Generating split_sequin.tsv
    writeQuin("split_sequin.tsv", stats, o);
    
    // Generating split_reads.tsv
    SWriteReads(Product::Genomics, "split_reads.tsv", stats, o);

    // Synthetic calibration
    StageTwoLadderCalibration(stats, GSplit::analyze, o.work, o, Standard::instance().gen.l3());

    // Generating ladder table
    writeLTable(o.work + "/split_ladder.tsv", "split_ladder_table.tsv", o);
    
    // Generating sequin abundance table
    writeSTable(o.work + "/split_sequin.tsv", "split_sequin_table.tsv", o, 6, 6, 6, "EXP_FREQ", "OBS_FREQ");

    // Generating split_summary.txt
    writeSummary("split_summary.txt", f1, f2, stats, o);
    
    if (o.report)
    {
        // Generating split_somatic.R
        writeSomatic("report_files/split_somatic.R", "split_sequin.tsv", o);

        // Generating HTML report
        GSplit::writeGReport("split_report.html", o);
    }
}
