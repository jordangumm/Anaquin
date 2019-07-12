#include "Meta/Meta.hpp"
#include "Meta/m_split.hpp"
#include "tools/report.hpp"
#include "data/resources.hpp"
#include "writers/r_writer.hpp"

using namespace Anaquin;

typedef MSplit::Stats Stats;
typedef MSplit::Options Options;

static void ladder(Stats &stats, const Options &o)
{
    const auto l1 = Standard::instance().meta.l1();
    const auto l3 = Standard::instance().meta.l3();

    for (const auto &x : stats.S1.K.seqs)
    {
        // Metagenomics ladder abundance
        if (MBin(x) == Bin::GR && stats.S1.R.d2u.meds.count(x) && l1->contains(x))
        {
            stats.l1.add(x, l1->input(x, o.mix), stats.S1.R.d2u.meds.at(x));
        }
        
        // Synthetic ladder
        if (MBin(x) == LD && l3->contains(x) && stats.S1.R.d2u.meds.count(x))
        {
            stats.l2.add(x, l3->input(x), stats.S1.R.d2u.meds.at(x));
        }
    }
}

Stats MSplit::analyze(const FileName &f1, const FileName &f2, const Options &o)
{
    assert((o.calibL == -1) || (o.calibL >= 0 && o.calibL <= 1.0));
    assert((o.calibS == -1) || (o.calibS >= 0 && o.calibS <= 1.0));

    o.info("Index: "     + o.index);
    o.info("Kmer: "      + toString(o.k));
    o.info("FASTQ: "     + toString(o.forceFQ));
    o.info("Skip: "      + toString(o.skipKM));
    o.info("Threads: "   + toString(o.thr));
    o.info("Threshold: " + toString(o.rule));
    o.info("Mixture: "   + std::string(o.mix == Mix_1 ? "A" : "B"));
    o.info("Calibration (sequins): " + (o.calibS == -1 ? "NA" : noTrailZero(std::to_string(o.calibS))));
    o.info("Calibration (ladder): "  + (o.calibL == -1 ? "NA" : noTrailZero(std::to_string(o.calibL))));

    Stats stats;

    // Kallisto on raw inputs before any calibration to "S1"
    SCombine(f1, stats.S1, o); SKallisto(stats.S1, f1, f2, o);

    // Calibrate "sequins" by percentage
    stats.C1 = SCalibrateP(GR, o.calibS, stats.S1, o, SCalibrateDefault("sequin"));
    
    // Make sure the working directory is unaffected
    auto o_ = cloneO(o);

    if (o.flip && !o.bam)
    {
        o_.flipBefore = true;
    }

    // Kallisto on calibrated "sequins" to "S2"
    SCombine(stats.C1.o1, stats.S2, o_); SKallisto(stats.S2, stats.C1.o1, stats.C1.o2, o_); removeD(o_.work);

    o.info("Checking ladder");
    ladder(stats, o);

    return stats;
}

static void writeSummary(const FileName &file, const FileName &f1, const FileName &f2, const Stats &stats, const Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "-------SUMMARY STATISTICS\n\n"
                        "-------REFERENCE FILES\n\n"
                        "       Reference index:   %1%\n"
                        "       Reference mixture: %2%\n"
                        "       Mixture type:      %3%\n\n"
                        "-------LIBRARY INFORMATION\n\n"
                        "       Version:           %4%\n"
                        "       Instrument ID:     %5%\n"
                        "       Run number:        %6%\n"
                        "       Flowcell ID:       %7%\n"
                        "       Lane:              %8%\n\n"
                        "-------USER-SUPPLIED FILES\n\n"
                        "       Input file (first):  %9%\n"
                        "       Input file (second): %10%\n\n"
                        "-------PARAMETERS\n\n"
                        "       K-mer length:       %11%\n"
                        "       Sequin Calibration: %12%\n"
                        "       Ladder Calibration: %13%\n"
                        "       Threshold:          %14%\n\n"
                        "-------PARTITION SUMMARY\n\n"
                        "       Sample reads: %15% (%16$.4f%%)\n"
                        "       Sequin reads: %17% (%18$.4f%%)\n"
                        "       Ladder reads: %19% (%20$.4f%%)\n"
                        "       Vector reads: %21% (%22$.4f%%)\n"
                        "       Info reads:   %23% (%24$.4f%%)\n"
                        "       Dilution:     %25%%%\n"
                        "       Total:        %26%\n\n"
                        "-------OUTPUT FASTQ FILES\n\n"
                        "       Sample reads path: %27%\n"
                        "       Sequin reads path: %28%\n"
                        "       Ladder reads path: %29%\n"
                        "       Vector reads path: %30%\n"
                        "       Info reads path:   %31%\n"
                        "       Calibrated sequin: %32%\n"
                        "       Calibrated ladder: %33%\n"
                        "       Merged:            %34%\n\n"
                        "-------CALIBRATION SUMMARY\n\n"
                        "       Sample (before): %35%\n"
                        "       Sequin (before): %36%\n"
                        "       Sequin target:   %37%\n"
                        "       Sequin scaling:  %38%\n"
                        "       Sequin (after):  %39%\n\n"
                        "       Ladder (before): %40%\n"
                        "       Ladder target:   %41%\n"
                        "       Ladder scaling:  %42%\n"
                        "       Ladder (after):  %43%\n\n"
                        "       Final calibrated library:  %44%\n"
                        "       Final calibrated dilution: %45%";
    
    #define C(x) (stats.S1.K.c1.count(x) ? stats.S1.K.c1.at(x) : 0)
    #define S2(x) toString(x,2)
    
    extern FASTQ __KFQInfo__;
    const auto f = __KFQInfo__.format();
    
    // Dilution after calibration
    const auto afterD = (stats.C1.aSeq + stats.C2.aSeq) / (float) (stats.C1.aSeq + stats.C2.aSeq + stats.C1.bSam);
    
    const auto s1 = (stats.C1.p == -1) ? "NA" : toString(stats.C1.p, 2);
    const auto s2 = (stats.C2.p == -1) ? "NA" : toString(stats.C2.p, 2);

    o.writer->write((boost::format(format) % o.index             // 1
                                           % Standard::instance().meta.l1()->src
                                           % (o.mix == Mix_1 ? "A" : "B")
                                           % SVersion(Standard::instance().meta, stats.S1.K) // 4
                                           % __KFQInfo__.inst(f) // 5
                                           % __KFQInfo__.run(f)  // 6
                                           % __KFQInfo__.flow(f) // 7
                                           % __KFQInfo__.lane(f) // 8
                                           % f1                  // 9
                                           % f2                  // 10
                                           % o.k                 // 11
                                           % ((o.calibS == -1) ? "NA" : noTrailZero(std::to_string(o.calibS)))
                                           % ((o.calibL == -1) ? "NA" : noTrailZero(std::to_string(o.calibL)))
                                           % o.rule              // 14
                                           % C(ES)               // 15
                                           % S2(100.0 * stats.S1.K.binP(ES))
                                           % C(GR)               // 17
                                           % S2(100.0 * stats.S1.K.binP(GR))
                                           % C(LD)               // 19
                                           % S2(100.0 * stats.S1.K.binP(LD))
                                           % C(VC)               // 21
                                           % S2(100.0 * stats.S1.K.binP(VC))
                                           % C(IF)               // 23
                                           % S2(100.0 * stats.S1.K.binP(IF))
                                           % toString(100.0 * stats.S1.dil(), 2)
                                           % stats.S1.K.total()         // 26
                                           % (o.work + "/meta_sample*") // 27
                                           % (o.work + "/meta_sequin*") // 28
                                           % (o.work + "/meta_ladder*") // 29
                                           % (o.work + "/meta_vector*") // 30
                                           % (o.work + "/meta_info*")   // 31
                                           % (o.work + "/meta_sequin_calibrated*") // 32
                                           % (o.work + "/meta_ladder_calibrated*") // 33
                                           % (o.work + "/meta_merged*") // 34
                                           % stats.C1.bSam              // 35
                                           % stats.C1.bSeq              // 36
                                           % stats.C1.tar               // 37
                                           % s1                         // 38
                                           % stats.C1.aSeq              // 39
                                           % stats.C2.bSeq              // 40
                                           % stats.C2.tar               // 41
                                           % s2                         // 42
                                           % stats.C2.aSeq              // 43
                                           % (stats.C1.aSeq + stats.C2.aSeq)
                                           % afterD                     // 45
                     ).str());
    o.writer->close();
}

static void writeQuins(const FileName &file, const SStats &stats, const Options &o)
{
    const auto form = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(form) % "Name"
                                         % "Mix"
                                         % "SD"
                                         % "Mean"
                                         % "Min"
                                         % "Q25"
                                         % "Med"
                                         % "Q75"
                                         % "Max"
                                         % "Count").str());
    
    #define _S1_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? "-" : toString(x.at(y)) : MISSING)
    #define _S2_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? "-" : toString(x.at(y), 2) : MISSING)

    const auto l1 = Standard::instance().meta.l1();
    assert(l1);
    
    for (const auto &x : stats.K.seqs)
    {
        const auto mix = (MBin(x) == LD || MBin(x) == IF || MBin(x) == VC) ? MISSING :
                         (l1->contains(x) ? toString(l1->input(x, o.mix), 4) : MISSING);
        const auto cn = stats.K.sqc.count(x) ? toString(stats.K.sqc.at(x)) : MISSING;
        
        o.writer->write((boost::format(form) % x
                                             % mix
                                             % _S1_(stats.R.d2u.sds, x)
                                             % _S1_(stats.R.d2u.mus, x)
                                             % _S1_(stats.R.d2u.mins, x)
                                             % _S2_(stats.R.d2u.q25,  x)
                                             % _S2_(stats.R.d2u.meds, x)
                                             % _S2_(stats.R.d2u.q75,  x)
                                             % _S1_(stats.R.d2u.maxs, x)
                                             % cn).str());
    }
    
    o.writer->close();
}

void MSplit::report(const FileName &f1, const FileName &f2, const Options &o)
{
    assert( Standard::instance().meta.l1());
    assert(!Standard::instance().meta.l2());
    assert( Standard::instance().meta.l3());
    
    auto stats = analyze(f1, f2, o);
    
    // Generating meta_sequin.tsv
    writeQuins("meta_sequin.tsv", stats.S1, o);

    if (o.calibS != -1)
    {
        // Generating meta_sequins_calibrated.tsv
        writeQuins("meta_sequins_calibrated.tsv", stats.S2, o);
    }

    const auto src = (o.calibS != -1) ? "meta_sequins_calibrated.tsv" : "meta_sequin.tsv";

    /*
     * Generating meta_abundance.R
     */
    
    o.generate("meta_abundance.R");
    o.writer->open("meta_abundance.R");
    o.writer->write(RWriter::createLinear(src, o.work, "Sequin Ladder", "Expected Abundance (Normalised, log2)", "Measured Abundance (Read count, log2)", "data$Mix", "data$Count"));
    o.writer->close();

    // Generating meta_reads.tsv
    SWriteReads(Product::Meta, "meta_reads.tsv", stats.S1, o);

    // Temporary for keeping the intermediate calibrated
    const auto o1 = cloneO(o);
    
    /*
     * Synthetic calibration from inputs ("stats.S1"). Read pre-calibration from "o.work" and write
     * results into "o1.work".
     */
    
    stats.S5 = SCalibSynthetic(stats.S1,
                   [&](const FileName &f1, const FileName &f2, const Options &o) {
                       return MSplit::analyze(f1, f2, o).S1; }, o.work, o1, Standard::instance().meta.l3());

    assert(stats.S3.K.total() == 0);
    
    if (!stats.S5.r1.empty())
    {
        // Check synthetic calibration to "S3"
        SCombine(stats.S5.r1, stats.S3, o1); SKallisto(stats.S3, stats.S5.r1, stats.S5.r2, o1);

        // Assign back to the sample reads for the next step
        stats.S3.K.c1[ES] = stats.S1.K.c1[ES];
        
        // Calibration by percentage for ladders, write results into output directory
        stats.C2 = SCalibrateP(LD, o.calibL, stats.S3, o, SCalibrateDefault("ladder", o1.work, o.work));
        
        // Kallisto on double calibrated "ladders" to "S4"
        const auto o2 = cloneO(o); SCombine(stats.C2.o1, stats.S4, o2); SKallisto(stats.S4, stats.C2.o1, stats.C2.o2, o2); removeD(o2.work);

        // This can only be done after "o2"
        removeD(o1.work);
        
        // Overwrite outputs generated by SCalibSynthetic()
        SWriteLadderPostCalib(Standard::instance().meta.l3(), stats.S4, o);
    }

    // Generating meta_summary.stats
    writeSummary("meta_summary.stats", f1, f2, stats, o);
    
    if (o.skipMerge)
    {
        o.info("Skipped merging");
    }
    else
    {
        o.info("Merging sample, sequin and ladder");

        if (o.writeBAM())
        {
            const auto seq = (o.calibS == -1) ? o.work + "/meta_sequin.bam" : o.work + "/meta_sequin_calibrated.bam";
            const auto lad = (o.calibL == -1) ? o.work + "/meta_ladder.bam" : o.work + "/meta_ladder_calibrated.bam";
            mergeBAM(std::vector<FileName> { o.work + "/meta_sample.bam", seq, lad }, o.work + "/meta_merged.bam");
        }
        else
        {
            const auto seq_1 = (o.calibS == -1) ? o.work + "/meta_sequin_1.fq.gz" : o.work + "/meta_sequin_calibrated_1.fq.gz";
            const auto seq_2 = (o.calibS == -1) ? o.work + "/meta_sample_2.fq.gz" : o.work + "/meta_sequin_calibrated_2.fq.gz";
            const auto lad_1 = (o.calibL == -1) ? o.work + "/meta_ladder_1.fq.gz" : o.work + "/meta_ladder_calibrated_1.fq.gz";
            const auto lad_2 = (o.calibL == -1) ? o.work + "/meta_ladder_2.fq.gz" : o.work + "/meta_ladder_calibrated_2.fq.gz";
            mergeFQ(std::vector<FileName> { o.work + "/meta_sample_1.fq.gz", seq_1, lad_1 },
                    std::vector<FileName> { o.work + "/meta_sample_2.fq.gz", seq_2, lad_2 },
                    o.work + "/meta_merged_1.fq.gz", o.work + "/meta_merged_2.fq.gz");
        }
    }
    
    // Generating HTML report
    if (o.report)
    {
        typedef Report::Options::KMOptions KOptions;
        Report::Options o2(o);
        o2.k = std::shared_ptr<KOptions>(new KOptions());
        Report::meta(o2);
    }
}
