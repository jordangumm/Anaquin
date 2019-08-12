#include "Meta/Meta.hpp"
#include "stats/stats.hpp"
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
    o.info("Index: "     + o.index);
    o.info("Kmer: "      + toString(o.k));
    o.info("FASTQ: "     + toString(o.forceFQ));
    o.info("Skip: "      + toString(o.skipKM));
    o.info("Threads: "   + toString(o.thr));
    o.info("Threshold: " + toString(o.rule));
    o.info("Mixture: "   + mixToStr(o.mix));
    o.info("Calibration (sequins): " + (o.oneS == NO_CALIBRATION ? MISSING : toString(o.oneS)));
    o.info("Calibration (ladder): "  + (!o.isLCalib() ? MISSING : toString(o.ladP())));

    Stats stats;

    // Kallisto on raw inputs before any calibration to "S1"
    SAlign(f1, stats.S1, o); SKallisto(stats.S1, f1, f2, o);

    // First-stage calibration on "sequins"
    stats.C1 = SCalibrateP(GR, o.oneS, stats.S1, o, SCalibrateDefault("sequin"));
    
    // Make sure the working directory is unaffected
    auto o_ = cloneO(o);

    if (o.flip && !o.bam)
    {
        o_.flipBefore = true;
    }

    // Kallisto on first-stage calibrated "sequins" to "S2"
    SAlign(stats.C1.o1, stats.S2, o_); SKallisto(stats.S2, stats.C1.o1, stats.C1.o2, o_); removeD(o_.work);

    ladder(stats, o);
    return stats;
}

static void writeSummary(const FileName &file, const FileName &f1, const FileName &f2, const FileName &tsv1, const FileName &tsv2, const Stats &stats, const Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto f = "SEQUIN REPORT:                 %1%\n\n"
                   "REFERENCE FILES\n"
                   "Reference index:               %2%\n"
                   "Reference mixture:             %3%\n"
                   "Mixture type:                  %4%\n\n"
                   "LIBRARY INFORMATION\n"
                   "Version:                       %5%\n"
                   "Instrument ID:                 %6%\n"
                   "Run number:                    %7%\n"
                   "Flowcell ID:                   %8%\n"
                   "Lane:                          %9%\n\n"
                   "USER-SUPPLIED FILES\n"
                   "Input file (first):            %10%\n"
                   "Input file (second):           %11%\n\n"
                   "ANAQUIN PARAMETERS\n"
                   "K-mer length:                  %12%\n"
                   "Threshold:                     %13%\n\n"
                   "PARTITION SUMMARY (BEFORE CALIBRATION)\n"
                   "Sample reads:                  %16% (%17%)\n"
                   "Sequin reads:                  %18% (%19%)\n"
                   "Ladder reads:                  %20% (%21%)\n"
                   "Vector reads:                  %22% (%23%)\n"
                   "Info reads:                    %24% (%25%)\n"
                   "Dilution:                      %26%%%\n"
                   "Total reads:                   %27%\n\n"
                   "OUTPUT FILES\n"
                   "Sample reads path:             %34%\n"
                   "Sequin reads path:             %35%\n"
                   "Ladder reads path:             %36%\n"
                   "Vector reads path:             %37%\n"
                   "Calibrated sequin:             %38%\n"
                   "Calibrated ladder:             %39%\n"
                   "Merged:                        %40%\n\n"
                   "CALIBRATION SUMMARY\n"
                   "Sequin calibration:            %14%\n"
                   "Sequin calibration factor:     %28%\n"
                   "Calibrated sequin reads:       %29%\n"
                   "Ladder calibration:            %15%\n"
                   "Ladder calibration factor:     %30%\n"
                   "Calibrated ladder reads:       %31%\n"
                   "Dilution after calibration:    %32%%%\n"
                   "Calibrated total reads:        %33%\n\n"
                   "LADDER - LIBRARY QUALITY\n"
                   "Slope:                         %41%\n"
                   "R2:                            %42%\n"
                   "Mean Ratio:                    %43%\n"
                   "Ladder table:                  %44%\n\n"
                   "SEQUIN QUANTIFICATION\n"
                   "Slope:                         %45%\n"
                   "R2:                            %46%\n"
                   "Sequin table:                  %47%";

    const auto tmp = tmpFile();
    
    // Only ladder sequins
    RGrep(tsv1, tmp, "NAME", "LD_"); const auto l1 = RLinear(tmp, "NAME", "UNIT", "READ").linear();
    
    // Only metagenomics sequins
    RGrep(tsv2, tmp, "NAME", "MG_"); const auto l2 = RLinear(tmp, "NAME", "MIX", "TPM").linear();

    #define C(x) (stats.S1.K.c1.count(x) ? stats.S1.K.c1.at(x) : 0)
    
    extern FASTQ __KFQInfo__;
    const auto fo = __KFQInfo__.format();
    
    /*
     * Total number of reads after calibration is sum of:
     *
     *    Number of calibrated sequin stage one reads
     *    Number of calibrated ladder stage one reads
     *    Number of calibrated ladder stage two reads
     */

    // Any calibration at all?
    const auto anyCalib = o.isSCalib() || o.isLCalib();
    
    // Calibrated total reads (only valid if there's calibration)
    const auto afterR = stats.C1.aSeq + stats.C2.aSeq + stats.S5.after;
    
    // Dilution after calibration
    const auto afterD = afterR / (float) (afterR + stats.C1.bSam);
    
    // Scaling factor for sequins
    const auto s1 = (stats.C1.p == NO_CALIBRATION) ? MISSING : toString(stats.C1.p, 2);

    // Missing if no ladder calibration
    auto t1 = MISSING; auto t2 = MISSING; auto t3 = NO_CALIBRATION;
    
    if (o.oneL != NO_CALIBRATION)
    {
        t1 = toString(o.oneL);
        t2 = toString(stats.C2.p);
        t3 = stats.C2.aSeq;
    }
    else if (o.secL != NO_CALIBRATION)
    {
        t1 = toString(o.secL);
        t2 = toString(stats.S5.meanS); // Different for each individual standard
        t3 = stats.S5.after;
    }
    
    o.writer->write((boost::format(f) % date()
                                      % o.index
                                      % Standard::instance().meta.l1()->src
                                      % mixToStr(o.mix)
                                      % SVersion(Standard::instance().meta, stats.S1.K)
                                      % __KFQInfo__.inst(fo) // 6
                                      % __KFQInfo__.run(fo)  // 7
                                      % __KFQInfo__.flow(fo) // 8
                                      % __KFQInfo__.lane(fo) // 9
                                      % f1                   // 10
                                      % f2                   // 11
                                      % o.k                  // 12
                                      % o.rule               // 13
                                      % (!o.isSCalib() ? MISSING : toString(o.oneS))
                                      % t1
                                      % C(ES)               // 16
                                      % S2(100.0 * stats.S1.K.binP(ES))
                                      % C(GR)               // 18
                                      % S2(100.0 * stats.S1.K.binP(GR))
                                      % C(LD)               // 20
                                      % S2(100.0 * stats.S1.K.binP(LD))
                                      % C(VC)               // 22
                                      % S2(100.0 * stats.S1.K.binP(VC))
                                      % C(IF)               // 24
                                      % S2(100.0 * stats.S1.K.binP(IF))
                                      % toString(100.0 * stats.S1.dil(), 2)
                                      % stats.S1.K.total()  // 27
                                      % s1                  // 28
                                      % (!o.isSCalib() ? MISSING : toString(stats.C1.aSeq)) // 29
                                      % t2                            // 30
                                      % (o.isLCalib() ? toString(t3) : MISSING) // 31
                                      % (anyCalib ? toString(100.0 * afterD) : MISSING) // 32
                                      % (anyCalib ? toString(afterR, 0) : MISSING)
                                      % (o.work + "/meta_sample*") // 34
                                      % (o.work + "/meta_sequin*") // 35
                                      % (o.work + "/meta_ladder*") // 36
                                      % (o.work + "/meta_vector*") // 37
                                      % (stats.C1.p == NO_CALIBRATION ? MISSING : o.work + "/meta_sequin_calibrated*") // 38
                                      % ((o.secL == NO_CALIBRATION && o.oneL == NO_CALIBRATION) ? MISSING : o.work + "/meta_ladder_calibrated*") // 39
                                      % (o.skipMerge ? MISSING : o.work + "/meta_merged*")
                                      % replaceNA(l1.m)  // 41
                                      % replaceNA(l1.R2) // 42
                                      % replaceNA(RLadTable(tsv1, tmp, "NAME"))              // 43
                                      % (o.work + "/meta_ladder_table.tsv") // 44
                                      % replaceNA(l2.m)  // 45
                                      % replaceNA(l2.R2) // 46
                                      % (o.work + "/meta_sequin_table.tsv") // 47
                     ).str());
    o.writer->close();
}

static void writeQuin(const FileName &file, const SStats &stats, const Options &o)
{
    const auto form = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%";
    
    #define _S1_(x,y) (x.count(y) && std::isnan(x.at(y)) ? toString(x.at(y)) : MISSING)
    #define _S2_(x,y) (x.count(y) && std::isnan(x.at(y)) ? toString(x.at(y), 2) : MISSING)
    
    const auto l1 = Standard::instance().meta.l1();
    assert(l1);

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(form) % "NAME"
                                         % "MIX"
                                         % "SD"
                                         % "MEAN"
                                         % "MIN"
                                         % "Q25"
                                         % "Q50"
                                         % "Q75"
                                         % "MAX"
                                         % "READ"
                                         % "TPM").str());

    // Unnormalized read counts
    auto count = [&](const SequinID &x)
    {
        return stats.K.sqc.count(x) ? stats.K.sqc.at(x) : std::numeric_limits<unsigned>::quiet_NaN();
    };

    // "per million" scaling factor
    auto scale = 0.0;
    
    /*
     * Sum up all the reads to derive scaling factor
     */
    
    for (const auto &x : stats.K.seqs)
    {
        if (!std::isnan(count(x)))
        {
            scale += count(x);
        }
    }

    // Divide by a million
    scale /= 1000000;

    for (const auto &x : stats.K.seqs)
    {
        const auto mix = (MBin(x) == LD || MBin(x) == IF || MBin(x) == VC) ? MISSING :
                         (l1->contains(x) ? toString(l1->input(x, o.mix), 6) : MISSING);

        // Raw count
        const auto cn = count(x);

        // Transcripts per million
        const auto tpm = !std::isnan(cn) && scale ? toString(cn / scale) : MISSING;
        
        o.writer->write((boost::format(form) % x
                                             % mix
                                             % _S1_(stats.R.d2u.sds, x)
                                             % _S1_(stats.R.d2u.mus, x)
                                             % _S1_(stats.R.d2u.mins, x)
                                             % _S2_(stats.R.d2u.q25,  x)
                                             % _S2_(stats.R.d2u.meds, x)
                                             % _S2_(stats.R.d2u.q75,  x)
                                             % _S1_(stats.R.d2u.maxs, x)
                                             % replaceNA(cn)
                                             % tpm).str());
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
    writeQuin("meta_sequin.tsv", stats.S1, o);

    if (o.isSCalib())
    {
        // Generating meta_sequin_calibrated.tsv
        writeQuin("meta_sequin_calibrated.tsv", stats.S2, o);
    }

    const auto src = (o.oneS != NO_CALIBRATION) ? "meta_sequin_calibrated.tsv" : "meta_sequin.tsv";

    // Generating meta_reads.tsv
    SWriteReads(Product::Meta, "meta_reads.tsv", stats.S1, o);

    // Temporary for keeping the intermediate calibrated
    const auto o1 = cloneO(o);
    
    /*
     * Stage two ladder calibration from Kallisto inputs ("stats.S1"). Read pre-calibration from "o.work" and write
     * results into "o1.work".
     */
    
    stats.S5 = StageTwoLadderCalibration(stats.S1,
                   [&](const FileName &f1, const FileName &f2, const Options &o) {
                       return MSplit::analyze(f1, f2, o).S1; }, o.work, o1, Standard::instance().meta.l3());

    assert(stats.S3.K.total() == 0);
    
    if (!stats.S5.r1.empty())
    {
        // Check stage two ladder calibration to "S3"
        SAlign(stats.S5.r1, stats.S3, o1); SKallisto(stats.S3, stats.S5.r1, stats.S5.r2, o1);

        // Assign back to the sample reads for the next step
        stats.S3.K.c1[ES] = stats.S1.K.c1[ES];
        
        /*
         * This is not a bug. Although it's first-stage, we apply it after the second-stage
         * in ladder calibration. Anaquin guarantees only one of those would applied.
         */
        
        // Calibration ladders, write results into output directory
        stats.C2 = SCalibrateP(LD, o.oneL, stats.S3, o, SCalibrateDefault("ladder", o1.work, o.work));
        
        // Kallisto on double calibrated (stage one and two, even only one would apply) "ladders" to "S4"
        const auto o2 = cloneO(o); SAlign(stats.C2.o1, stats.S4, o2); SKallisto(stats.S4, stats.C2.o1, stats.C2.o2, o2); removeD(o2.work);

        // This can only be done after "o2"
        removeD(o1.work);
        
        // Overwrite outputs generated by StageTwoLadderCalibration()
        SWriteLadderPostCalib(Standard::instance().meta.l3(), stats.S4, o);
    }
    
    const auto tsv1 = o.isLCalib() ? o.work + "/meta_ladder_calibrated.tsv" : o.work + "/meta_ladder.tsv";
    const auto tsv2 = o.isSCalib() ? o.work + "/meta_sequin_calibrated.tsv" : o.work + "/meta_sequin.tsv";
    
    // Generating ladder table
    writeLTable(tsv1, "meta_ladder_table.tsv", o);
    
    // Generating sequin abundance table
    writeSTable(tsv2, "meta_sequin_table.tsv", o, 6, 6, 6, "MIX", "TPM");

    // Generating meta_summary.txt
    writeSummary("meta_summary.txt", f1, f2, tsv1, tsv2, stats, o);
    
    if (o.skipMerge)
    {
        o.info("Skipped merging");
    }
    else
    {
        o.info("Merging sample, sequin and ladder reads");

        if (o.writeBAM())
        {
            const auto seq = (!o.isSCalib()) ? o.work + "/meta_sequin.bam" : o.work + "/meta_sequin_calibrated.bam";
            const auto lad = (!o.isLCalib()) ? o.work + "/meta_ladder.bam" : o.work + "/meta_ladder_calibrated.bam";
            mergeBAM(std::vector<FileName> { o.work + "/meta_sample.bam", seq, lad }, o.work + "/meta_merged.bam");
        }
        else
        {
            const auto seq_1 = (!o.isSCalib()) ? o.work + "/meta_sequin_1.fq.gz" : o.work + "/meta_sequin_calibrated_1.fq.gz";
            const auto seq_2 = (!o.isSCalib()) ? o.work + "/meta_sample_2.fq.gz" : o.work + "/meta_sequin_calibrated_2.fq.gz";
            const auto lad_1 = (!o.isLCalib()) ? o.work + "/meta_ladder_1.fq.gz" : o.work + "/meta_ladder_calibrated_1.fq.gz";
            const auto lad_2 = (!o.isLCalib()) ? o.work + "/meta_ladder_2.fq.gz" : o.work + "/meta_ladder_calibrated_2.fq.gz";
            mergeFQ(std::vector<FileName> { o.work + "/meta_sample_1.fq.gz", seq_1, lad_1 },
                    std::vector<FileName> { o.work + "/meta_sample_2.fq.gz", seq_2, lad_2 },
                    o.work + "/meta_merged_1.fq.gz", o.work + "/meta_merged_2.fq.gz");
        }
    }
    
    // Generating HTML report
    if (o.report)
    {
        /*
         * Generating meta_abundance.R
         */
        
        o.generate("meta_abundance.R");
        o.writer->open("report_files/meta_abundance.R");
        o.writer->write(RWriter::createLinear(src, o.work, "Sequin Ladder", "Expected Abundance (Normalised, log2)", "Measured Abundance (TPM, log2)", "data$MIX", "data$TPM"));
        o.writer->close();

        typedef Report::Options::KMOptions KOptions;
        Report::Options o2(o);
        o2.k = std::shared_ptr<KOptions>(new KOptions());
        Report::meta(o2);
    }
}
