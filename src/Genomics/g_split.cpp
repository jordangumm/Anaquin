#include <libgen.h>
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
                
                stats.VR[seq] = hasR ? SS::med(toVector(stats.K.uniqs.at(R))) : NAN;
                stats.VV[seq] = hasV ? SS::med(toVector(stats.K.uniqs.at(V))) : NAN;
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
                const auto rn = stats.K.uniqs.count(R) ? SS::med(toVector(stats.K.uniqs.at(R))) : 0;
                
                // Median k-mer for variant sequin
                const auto vn = stats.K.uniqs.count(V) ? SS::med(toVector(stats.K.uniqs.at(V))) : 0;
                
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
    SCombine(f1, stats, o, nullptr, true); SKallisto(stats, f1, f2, o);

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

    const auto format = "-------SUMMARY STATISTICS\n\n"
                        "-------REFERENCE FILES\n\n"
                        "       Reference index: %1%\n\n"
                        "-------LIBRARY INFORMATION\n\n"
                        "       Version:       %43%\n"
                        "       Instrument ID: %44%\n"
                        "       Run number:    %45%\n"
                        "       Flowcell ID:   %46%\n"
                        "       Lane:          %47%\n\n"
                        "-------USER-SUPPLIED FILES\n\n"
                        "       Input file (first):  %2%\n"
                        "       Input file (second): %3%\n\n"
                        "-------PARAMETERS\n\n"
                        "       K-mer length: %4%\n"
                        "       Threshold:    %5%\n\n"
                        "-------PARTITION SUMMARY\n\n"
                        "       Human reads:             %6%  (%7$.4f%%)\n"
                        "       Genome reads:            %8% (%9$.4f%%)\n"
                        "       Ladder reads:            %10% (%11$.4f%%)\n"
                        "       Vector reads:            %12% (%13$.4f%%)\n"
                        "       Structural reads:        %14% (%15$.4f%%)\n"
                        "       Immune reads:            %16% (%17$.4f%%)\n"
                        "       HLA reads:               %18% (%19$.4f%%)\n"
                        "       Information reads:       %20% (%21$.4f%%)\n"
                        "       Mitochondria reads:      %22% (%23$.4f%%)\n"
                        "       Dilution:                %24%%%\n"
                        "       Total:                   %25%\n\n"
                        "-------OUTPUT FASTQ FILES\n\n"
                        "       Human reads path:        %26%\n"
                        "       Genome reads path:       %27%\n"
                        "       Vector reads path:       %28%\n"
                        "       Ladder reads path:       %29%\n"
                        "       Structural reads path:   %30%\n"
                        "       Immune reads path:       %31%\n"
                        "       HLA reads path:          %32%\n"
                        "       Information reads path:  %33%\n"
                        "       Mitochondria reads path: %34%\n\n"
                        "-------Somatic Allele Frequency Ladder (log2 scale)\n\n"
                        "       Slope (AF):       %35%\n"
                        "       R2 (AF):          %36%\n"
                        "       F-statistic: (AF) %37%\n"
                        "       P-value (AF):     %38%\n\n"
                        "-------Synthetic Ladder (log2 scale)\n\n"
                        "       Slope (SL):       %39%\n"
                        "       R2 (SL):          %40%\n"
                        "       F-statistic (SL): %41%\n"
                        "       P-value (SL):     %42%\n";

    LinearModel af, ld;
    try { af = stats.af.linear(true, true); } catch(...) {}
    try { ld = stats.ld.linear(true, true); } catch(...) {}

    const auto &r = Standard::instance().gen;
    
    #define C(x) (stats.K.c1.count(x) ? stats.K.c1.at(x) : 0)
    #define P(x) (100.0 * stats.K.binP(x))

    const auto gn = C(GR) + C(SO) + C(MS) + C(HP);
    const auto gp = P(GR) + P(SO) + P(MS) + P(HP);
    assert(gp >= 0 && gp <= 100);
    
    extern FASTQ __KFQInfo__;
    const auto f = __KFQInfo__.format();

    o.writer->write((boost::format(format) % o.index         // 1
                                           % f1              // 2
                                           % f2              // 3
                                           % o.k             // 4
                                           % o.rule          // 5
                                           % S0(C(ES))       // 6
                                           % S2(P(ES))       // 7
                                           % S0(gn)          // 8
                                           % S2(gp)          // 9
                                           % S0(C(LD))       // 10
                                           % S2(P(LD))       // 11
                                           % S0(C(VC))       // 12
                                           % S2(P(VC))       // 13
                                           % S0(C(SV))       // 14
                                           % S2(P(SV))       // 15
                                           % S0(C(IM))       // 16
                                           % S2(P(IM))       // 17
                                           % S0(C(HL))       // 18
                                           % S2(P(HL))       // 19
                                           % S0(C(IF))       // 20
                                           % S2(P(IF))       // 21
                                           % S0(C(MT))       // 22
                                           % S2(P(MT))       // 23
                                           % (100.0 * stats.dil())
                                           % stats.K.total()
                                           % (o.work + "/split_sample*") // 26
                                           % (o.work + "/split_sequin*") // 27
                                           % (o.work + "/split_vector*") // 28
                                           % (o.work + "/split_ladder*") // 29
                                           % (o.work + "/split_sv*")     // 30
                                           % (o.work + "/split_immune*") // 31
                                           % (o.work + "/split_hla*")    // 32
                                           % (o.work + "/split_info*")   // 33
                                           % (o.work + "/split_mito*")   // 34
                                           % MISS_STRING(af.m)           // 35
                                           % MISS_STRING(af.r)           // 36
                                           % MISS_STRING(af.F)           // 37
                                           % MISS_STRING(af.p)           // 38
                                           % MISS_STRING(ld.m)           // 39
                                           % MISS_STRING(ld.r)           // 40
                                           % MISS_STRING(ld.F)           // 41
                                           % MISS_STRING(ld.p)           // 42
                                           % SVersion(Standard::instance().gen, stats.K) // 43
                                           % __KFQInfo__.inst(f)                         // 44
                                           % __KFQInfo__.run(f)                          // 45
                                           % __KFQInfo__.flow(f)                         // 46
                                           % __KFQInfo__.lane(f)                         // 47
                     ).str());
    o.writer->close();
}

void GSplit::writeQuin(const FileName &file, const Stats &stats, const SOptions &o)
{
    const auto &r = Standard::instance().gen;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%\t%14%\t%15%\t%16%%17%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Label"
                                           % "Chrom"
                                           % "Position"
                                           % "ExpFreq"
                                           % "ObsFreq"
                                           % "RefCount"
                                           % "VarCount"
                                           % "Min"
                                           % "Q25"
                                           % "Med"
                                           % "Q75"
                                           % "Max"
                                           % "Read"
                                           % "Genotype"
                                           % "Type"
                                           % r.a1()->strForKeys()).str());
    
    std::map<SequinID, Base> r1;
    ParserFA::size(o.index, r1);

    #define _S1_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? "-" : toString(x.at(y)) : MISSING)
    #define _S2_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? "-" : toString(x.at(y), 2) : MISSING)

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
                                                   % removeNaN(exp)
                                                   % removeNaN(obs)
                                                   % removeNaN(R)
                                                   % removeNaN(V)
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
    
    // Generating split_summary.stats
    writeSummary("split_summary.stats", f1, f2, stats, o);

    // Generating split_sequin.tsv
    writeQuin("split_sequin.tsv", stats, o);
    
    // Generating split_somatic.R
    writeSomatic("split_somatic.R", "split_sequin.tsv", o);

    // Generating split_reads.tsv
    SWriteReads(Product::Genomics, "split_reads.tsv", stats, o);

    // Synthetic calibration
    SCalibSynthetic(stats, GSplit::analyze, o.work, o, Standard::instance().gen.l3());

    // Generating HTML report
    GSplit::writeGReport("split_report.html", o);
}
