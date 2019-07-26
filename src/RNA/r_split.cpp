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
    o.info("Index: "     + o.index);
    o.info("Kmer: "      + toString(o.k));
    o.info("FASTQ: "     + toString(o.forceFQ));
    o.info("Skip: "      + toString(o.skipKM));
    o.info("Threads: "   + toString(o.thr));
    o.info("Threshold: " + toString(o.rule));
    o.info("Mixture: "   + mixToStr(o.mix));
    o.info("Calibration (sequins): " + (!o.isSCalib() ? MISSING : toString(o.oneS)));
    
    Stats stats;
    
    // Kallisto on raw inputs before any calibration to "S1"
    SAlign(f1, stats.S1, o); SKallisto(stats.S1, f1, f2, o);
    
    // First stage calibration on sequins
    stats.C1 = SCalibrateP(GR, o.oneS, stats.S1, o, SCalibrateDefault("sequin"));
    
    // Make sure the working directory is unaffected
    auto o_ = cloneO(o);
    
    // Kallisto on calibrated "sequins" to "S2"
    SAlign(stats.C1.o1, stats.S2, o_); SKallisto(stats.S2, stats.C1.o1, stats.C1.o2, o_); removeD(o_.work);
    
    return stats;
}

static void writeSummary(const FileName &file, const FileName &f1, const FileName &f2, const Stats &stats, const FileName &tsv, const Options &o)
{
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
                   "Sample reads:                  %14% (%15%)\n"
                   "Sequin reads:                  %16% (%17%)\n"
                   "Dilution:                      %18%%%\n"
                   "Total reads:                   %19%\n\n"
                   "OUTPUT FILES\n"
                   "Sample reads path:             %20%\n"
                   "Sequin reads path:             %21%\n"
                   "Calibrated sequin:             %22%\n"
                   "Merged:                        %23%\n\n"
                   "CALIBRATION SUMMARY\n"
                   "Sequin calibration:            %24%\n"
                   "Sequin calibration factor:     %25%\n"
                   "Calibrated sequin reads:       %26%\n"
                   "Dilution after calibration:    %27%%%\n"
                   "Calibrated total reads:        %28%\n\n"
                   "SEQUIN QUANTIFICATION\n"
                   "Slope:                         %29%\n"
                   "R2:                            %30%\n"
                   "Isoform expression table:      %31%\n\n"
                   "GENE QUANTIFICATION\n"
                   "Slope:                         %32%\n"
                   "R2:                            %33%\n"
                   "Gene expresssion table:        %34%";

    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();

    // Sequins at the isoform level
    RGrep(tsv, tmp1, "NAME", "R"); const auto l1 = RLinear(tmp1, "NAME", "MIX", "TPM").linear();

    // Generating sequin abundance table
    writeSTable(tmp1, "rna_sequin_isoform_table.tsv", o, 6, 6, 6, "MIX", "TPM");
    
    RFilterC(tmp1, tmp2, std::set<Label> { "GENE", "MIX", "TPM" }, true);
    
    // Sequins at the gene level
    RAggregateSum(tmp2, tmp1, "GENE", Imputation::ToZero); const auto l2 = RLinear(tmp1, "GENE", "MIX", "TPM").linear();

    // Generating sequin abundance table
    writeSTable(tmp1, "rna_sequin_gene_table.tsv", o, 6, 6, 6, "MIX", "TPM");
    
    #define C(x) (stats.S1.K.c1.count(x) ? stats.S1.K.c1.at(x) : 0)
    
    extern FASTQ __KFQInfo__;
    const auto fo = __KFQInfo__.format();

    // Any calibration at all?
    const auto anyCalib = o.isSCalib();

    // Calibrated total reads (only valid if there's calibration)
    const auto afterR = stats.C1.aSeq;

    // Dilution after calibration
    const auto afterD = (stats.C1.aSeq) / (float) (stats.C1.aSeq + stats.C1.bSam);
    
    const auto s1 = !o.isSCalib() ? "NA" : toString(stats.C1.p, 2);
    
    o.generate(file);
    o.writer->open(file);    
    o.writer->write((boost::format(f) % date()
                                      % o.index             // 2
                                      % Standard::instance().rna.l1()->src
                                      % mixToStr(o.mix)
                                      % SVersion(Standard::instance().rna, stats.S1.K)
                                      % __KFQInfo__.inst(fo) // 6
                                      % __KFQInfo__.run(fo)  // 7
                                      % __KFQInfo__.flow(fo) // 8
                                      % __KFQInfo__.lane(fo) // 9
                                      % f1                   // 10
                                      % f2                   // 11
                                      % o.k                  // 12
                                      % o.rule               // 13
                                      % C(ES)                // 14
                                      % S2(100.0 * stats.S1.K.binP(ES))
                                      % C(GR)                // 16
                                      % S2(100.0 * stats.S1.K.binP(GR))
                                      % toString(100.0 * stats.S1.dil(), 2)
                                      % stats.S1.K.total()  // 19
                                      % (o.work + "/rna_sample*") // 20
                                      % (o.work + "/rna_sequin*") // 21
                                      % (stats.C1.p == NO_CALIBRATION ? MISSING : o.work + "/rna_sequin_calibrated*") // 22
                                      % (o.skipMerge ? MISSING : o.work + "/rna_merged*")
                                      % (!o.isSCalib() ? MISSING : toString(o.oneS)) // 24
                                      % s1                  // 25
                                      % (!o.isSCalib()  ? MISSING : toString(stats.C1.aSeq))
                                      % (anyCalib ? toString(100.0 * afterD) : MISSING) // 27
                                      % (anyCalib ? toString(afterR, 0) : MISSING) // 28
                                      % replaceNA(l1.m)  // 29
                                      % replaceNA(l1.R2) // 30
                                      % (o.work + "/rna_sequin_isoform_table.tsv") // 31
                                      % replaceNA(l2.m)  // 32
                                      % replaceNA(l2.R2) // 33
                                      % (o.work + "/rna_sequin_gene_table.tsv") // 34
                     ).str());
    o.writer->close();
}

static void writeQuin(const FileName &file, const SStats &stats, const Options &o)
{
    const auto form = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%\t%13%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(form) % "NAME"
                                         % "GENE"
                                         % "LENGTH"
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
    
    #define _S1_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? MISSING : toString(x.at(y)) : MISSING)
    #define _S2_(x,y) (x.count(y) ? std::isnan(x.at(y)) ? MISSING : toString(x.at(y), 2) : MISSING)
    
    const auto l1 = Standard::instance().rna.l1();
    const auto l3 = Standard::instance().rna.l3();
    
    assert(l1);
    assert(l3);

    // "per million" scaling factor
    auto scale = 0.0;

    /*
     * Sum up all the reads to derive scaling factor
     */
    
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
    o.writer->open("report_files/rna_gene.R");
    o.writer->write(RWriter::createGene(src, o.work,
                                        "Gene Expression",
                                        "Expected Abundance (Normalised, log2)",
                                        "Measured Abundance (TPM, log2)",
                                        "data$MIX",
                                        "data$TPM"));
    o.writer->close();
}

void RSplit::report(const FileName &f1, const FileName &f2, const Options &o)
{
    assert(Standard::instance().rna.l1());
    assert(Standard::instance().rna.l2());
    
    auto stats = analyze(f1, f2, o);
    
    // Generating rna_sequin.tsv
    writeQuin("rna_sequin.tsv", stats.S1, o);
    
    if (o.oneS != -1)
    {
        // Generating rna_sequin_calibrated.tsv
        writeQuin("rna_sequin_calibrated.tsv", stats.S2, o);
    }

    const auto src = (o.oneS != -1) ? "rna_sequin_calibrated.tsv" : "rna_sequin.tsv";

    // Generating rna_reads.tsv
    SWriteReads(Product::Meta, "rna_reads.tsv", stats.S1, o);

    const auto tsv = o.isSCalib() ? o.work + "/rna_sequin_calibrated.tsv" : o.work + "/rna_sequin.tsv";
    
    // Generating rna_summary.txt
    writeSummary("rna_summary.txt", f1, f2, stats, tsv, o);
    
    if (o.skipMerge)
    {
        o.info("Skipped merging");
    }
    else
    {
        o.info("Merging sample and sequin");
        
        if (o.writeBAM())
        {
            const auto seq = !o.isSCalib() ? o.work + "/rna_sequin.bam" : o.work + "/rna_sequin_calibrated.bam";
            mergeBAM(std::vector<FileName> { o.work + "/rna_sample.bam", seq }, o.work + "/rna_merged.bam");
        }
        else
        {
            const auto seq_1 = !o.isSCalib() ? o.work + "/rna_sequin_1.fq.gz" : o.work + "/rna_sequin_calibrated_1.fq.gz";
            const auto seq_2 = !o.isSCalib() ? o.work + "/rna_sample_2.fq.gz" : o.work + "/rna_sequin_calibrated_2.fq.gz";
            mergeFQ(std::vector<FileName> { o.work + "/rna_sample_1.fq.gz", seq_1  },
                    std::vector<FileName> { o.work + "/rna_sample_2.fq.gz", seq_2  },
                    o.work + "/rna_merged_1.fq.gz", o.work + "/rna_merged_2.fq.gz");
        }
    }
    
    if (o.report)
    {
        /*
         * Generating rna_isoform.R
         */
        
        o.generate("rna_isoform.R");
        o.writer->open("report_files/rna_isoform.R");
        o.writer->write(RWriter::createLinear(src, o.work, "Isoform Expression", "Expected Abundance (Normalised, log2)", "Measured Abundance (TPM, log2)", "data$MIX", "data$TPM"));
        o.writer->close();

        // Generating rna_gene.R
        createGene(src, o);

        typedef Report::Options::KMOptions KOptions;
        Report::Options o2(o);
        o2.k = std::shared_ptr<KOptions>(new KOptions());
        Report::rna(o2);
    }
}
