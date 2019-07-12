#include "RNA/RNA.hpp"
#include "RNA/r_split.hpp"
#include "tools/report.hpp"
#include "data/resources.hpp"
#include "writers/r_writer.hpp"

using namespace Anaquin;

typedef RSplit::Stats Stats;
typedef RSplit::Options Options;

Stats RSplit::analyze(const FileName &f1, const FileName &f2, const Options &o)
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
    o.info("Calibration (sequins): " + (o.calibS == -1 ? "-" : noTrailZero(std::to_string(o.calibS))));
    
    Stats stats;
    
    // Kallisto on raw inputs before any calibration to "S1"
    SCombine(f1, stats.S1, o); SKallisto(stats.S1, f1, f2, o);
    
    // Calibrate "sequins" by percentage
    stats.C1 = SCalibrateP(GR, o.calibS, stats.S1, o, SCalibrateDefault("sequin"));
    
    // Make sure the working directory is unaffected
    auto o_ = cloneO(o);
    
    // Kallisto on calibrated "sequins" to "S2"
    SCombine(stats.C1.o1, stats.S2, o_); SKallisto(stats.S2, stats.C1.o1, stats.C1.o2, o_); removeD(o_.work);
    
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
                        "       Threshold:          %13%\n\n"
                        "-------PARTITION SUMMARY\n\n"
                        "       Sample reads: %14% (%15$.4f%%)\n"
                        "       Sequin reads: %16% (%17$.4f%%)\n"
                        "       Vector reads: %18% (%19$.4f%%)\n"
                        "       Info reads:   %20% (%21$.4f%%)\n"
                        "       Dilution:     %22%%%\n"
                        "       Total:        %23%\n\n"
                        "-------OUTPUT FASTQ FILES\n\n"
                        "       Sample reads path: %24%\n"
                        "       Sequin reads path: %25%\n"
                        "       Vector reads path: %26%\n"
                        "       Info reads path:   %27%\n"
                        "       Calibrated sequin: %28%\n"
                        "       Merged:            %29%\n\n"
                        "-------CALIBRATION SUMMARY\n\n"
                        "       Sample (before): %30%\n"
                        "       Sequin (before): %31%\n"
                        "       Sequin target:   %32%\n"
                        "       Sequin scaling:  %33%\n"
                        "       Sequin (after):  %34%\n\n"
                        "       Final calibrated library:  %35%\n"
                        "       Final calibrated dilution: %36%";
    
    #define C(x) (stats.S1.K.c1.count(x) ? stats.S1.K.c1.at(x) : 0)
    #define S2(x) toString(x,2)
    
    extern FASTQ __KFQInfo__;
    const auto f = __KFQInfo__.format();
    
    // Dilution after calibration
    const auto afterD = (stats.C1.aSeq) / (float) (stats.C1.aSeq + stats.C1.bSam);
    
    const auto s1 = (stats.C1.p == -1) ? "NA" : toString(stats.C1.p, 2);
    
    o.writer->write((boost::format(format) % o.index             // 1
                                           % Standard::instance().rna.l1()->src
                                           % (o.mix == Mix_1 ? "A" : "B")
                                           % "-"                 // 4
                                           % __KFQInfo__.inst(f) // 5
                                           % __KFQInfo__.run(f)  // 6
                                           % __KFQInfo__.flow(f) // 7
                                           % __KFQInfo__.lane(f) // 8
                                           % f1                  // 9
                                           % f2                  // 10
                                           % o.k                 // 11
                                           % ((o.calibS == -1) ? "NA" : noTrailZero(std::to_string(o.calibS)))
                                           % o.rule              // 13
                                           % C(ES)               // 14
                                           % S2(100.0 * stats.S1.K.binP(ES))
                                           % C(GR)               // 16
                                           % S2(100.0 * stats.S1.K.binP(GR))
                                           % C(VC)               // 18
                                           % S2(100.0 * stats.S1.K.binP(VC))
                                           % C(IF)               // 20
                                           % S2(100.0 * stats.S1.K.binP(IF))
                                           % toString(100.0 * stats.S1.dil(), 2)
                                           % stats.S1.K.total()        // 23
                                           % (o.work + "/rna_sample*") // 24
                                           % (o.work + "/rna_sequin*") // 25
                                           % (o.work + "/rna_vector*") // 26
                                           % (o.work + "/rna_info*")   // 27
                                           % (o.work + "/rna_sequin_calibrated*") // 28
                                           % (o.work + "/rna_merged*")  // 29
                                           % stats.C1.bSam              // 30
                                           % stats.C1.bSeq              // 31
                                           % stats.C1.tar               // 32
                                           % s1                         // 33
                                           % stats.C1.aSeq              // 34
                                           % (stats.C1.bSam + stats.C1.aSeq)
                                           % afterD                     // 36
                     ).str());
    o.writer->close();
}

static void writeQuins(const FileName &file, const SStats &stats, const Options &o)
{
    const auto form = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(form) % "Name"
                                         % "Gene"
                                         % "Length"
                                         % "Mix"
                                         % "SD"
                                         % "Mean"
                                         % "Min"
                                         % "Q25"
                                         % "Med"
                                         % "Q75"
                                         % "Max"
                                         % "Count"
                                         % "TPM").str());
    
    #define _S1_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? "-" : toString(x.at(y)) : MISSING)
    #define _S2_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? "-" : toString(x.at(y), 2) : MISSING)
    
    const auto l1 = Standard::instance().rna.l1();
    const auto l3 = Standard::instance().rna.l3();
    
    assert(l1);
    assert(l3);

    // "per million" scaling factor
    auto scale = 0.0;
    
    for (const auto &x : stats.K.seqs)
    {
        assert(l3->contains(x) && l3->input(x));
        
        // Raw read count
        const auto read = stats.K.sqc.count(x) ? stats.K.sqc.at(x) : 0;
        
        // Reads per kilobase (RPK)
        const auto rpk = read / l3->input(x);

        scale += rpk;
    }
    
    // Divide by a million
    scale /= 1000000;
    
    for (const auto &x : stats.K.seqs)
    {
        const auto mix = (RBin(x) == LD || RBin(x) == IF || RBin(x) == VC) ? MISSING :
                         (l1->contains(x) ? std::to_string(l1->input(x, o.mix)) : MISSING);

        const auto r1 = stats.K.sqc.count(x) ? stats.K.sqc.at(x) : 0;
        const auto r2 = stats.K.sqc.count(x) ? toString(stats.K.sqc.at(x)) : MISSING;
        
        // Reads per kilobase (RPK)
        const auto rpk = r1 / l3->input(x);
        
        // Transcripts per million
        const auto tpm = scale ? toString(rpk / scale) : MISSING;

        o.writer->write((boost::format(form) % x
                                             % noLast(x, "_")
                                             % l3->input(x)
                                             % mix
                                             % _S1_(stats.R.d2u.sds, x)
                                             % _S1_(stats.R.d2u.mus, x)
                                             % _S1_(stats.R.d2u.mins, x)
                                             % _S2_(stats.R.d2u.q25,  x)
                                             % _S2_(stats.R.d2u.meds, x)
                                             % _S2_(stats.R.d2u.q75,  x)
                                             % _S1_(stats.R.d2u.maxs, x)
                                             % r2
                                             % (r2 == MISSING ? MISSING : tpm)).str());
    }
    
    o.writer->close();
}

static void createGene(const FileName &src, const RSplit::Options &o)
{
    o.generate("rna_gene.R");
    o.writer->open("rna_gene.R");
    o.writer->write(RWriter::createGene(src, o.work,
                                        "Gene Expression",
                                        "Expected Abundance (Normalised, log2)",
                                        "Measured Abundance (TPM, log2)",
                                        "data$Mix",
                                        "data$TPM"));
    o.writer->close();
}

void RSplit::report(const FileName &f1, const FileName &f2, const Options &o)
{
    assert(Standard::instance().rna.l1());
    assert(Standard::instance().rna.l2());
    
    auto stats = analyze(f1, f2, o);
    
    // Generating rna_sequin.tsv
    writeQuins("rna_sequin.tsv", stats.S1, o);
    
    if (o.calibS != -1)
    {
        // Generating rna_sequins_calibrated.tsv
        writeQuins("rna_sequins_calibrated.tsv", stats.S2, o);
    }

    const auto src = (o.calibS != -1) ? "rna_sequins_calibrated.tsv" : "rna_sequin.tsv";

    /*
     * Generating rna_isoform.R
     */
    
    o.generate("rna_isoform.R");
    o.writer->open("rna_isoform.R");
    o.writer->write(RWriter::createLinear(src, o.work, "Isoform Expression", "Expected Abundance (Normalised, log2)", "Measured Abundance (TPM, log2)", "data$Mix", "data$TPM"));
    o.writer->close();
    
    // Generating rna_gene.R
    createGene(src, o);

    // Generating rna_reads.tsv
    SWriteReads(Product::Meta, "rna_reads.tsv", stats.S1, o);
    
    // Generating rna_summary.stats
    writeSummary("rna_summary.stats", f1, f2, stats, o);
    
    if (o.skipMerge)
    {
        o.info("Skipped merging");
    }
    else
    {
        o.info("Merging sample and sequin");
        
        if (o.writeBAM())
        {
            const auto seq = (o.calibS == -1) ? o.work + "/rna_sequin.bam" : o.work + "/rna_sequin_calibrated.bam";
            mergeBAM(std::vector<FileName> { o.work + "/rna_sample.bam", seq }, o.work + "/rna_merged.bam");
        }
        else
        {
            const auto seq_1 = (o.calibS == -1) ? o.work + "/rna_sequin_1.fq.gz" : o.work + "/rna_sequin_calibrated_1.fq.gz";
            const auto seq_2 = (o.calibS == -1) ? o.work + "/rna_sample_2.fq.gz" : o.work + "/rna_sequin_calibrated_2.fq.gz";
            mergeFQ(std::vector<FileName> { o.work + "/rna_sample_1.fq.gz", seq_1  },
                    std::vector<FileName> { o.work + "/rna_sample_2.fq.gz", seq_2  },
                    o.work + "/rna_merged_1.fq.gz", o.work + "/rna_merged_2.fq.gz");
        }
    }
    
    // Generating HTML report
    if (o.report)
    {
        typedef Report::Options::KMOptions KOptions;
        Report::Options o2(o);
        o2.k = std::shared_ptr<KOptions>(new KOptions());
        Report::rna(o2);
    }
}
