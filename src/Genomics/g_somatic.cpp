#include <cmath>
#include "tools/report.hpp"
#include "Genomics/g_somatic.hpp"

using namespace Anaquin;

template <typename T, typename O> void writeQualFilter(const T &, const O &o)
{
    o.generate("somatic_qualFilter.R");
    o.writer->open("report_files/somatic_qualFilter.R");
    o.writer->write((boost::format(PlotQualFilter()) % date()
                                                     % o.cmd
                                                     % o.work
                                                     % "somatic_sequin.tsv").str());
    o.writer->close();
}

template <typename O> void writeAllele(const O &o)
{
    o.generate("somatic_ladder.R");
    o.writer->open("report_files/somatic_ladder.R");
    o.writer->write((boost::format(PlotAllele()) % date()
                                                 % o.cmd
                                                 % o.work
                                                 % "somatic_sequin.tsv").str());
    o.writer->close();
}

static void writeTable(const FileName &file, const GSomatic::Stats &stats, const GVariant::Options &o)
{
    const auto tmp1  = tmpFile();
    const auto tmp2  = tmpFile();
    const auto src   = o.work + "/somatic_sequin.tsv";
    const auto src1  = tmpFile(); // Filtered source including all detected sequins
    const auto src2  = tmpFile(); // Filtered source but only CM sequins
    const auto src1_ = tmpFile(); // src1 but only the interesting columns
    const auto src2_ = tmpFile(); // src2 but only the interesting columns
    const auto tmpM  = tmpFile(); // Temporary for means
    const auto tmpS  = tmpFile(); // Temporary for standard deviation

    RGrep(src,  tmp1, "LABEL", "SV", false); // Remove all sample variants not giving allele frequency
    RGrep(tmp1, tmp2, "LABEL", "FP", false); // Remove all false positives not giving allele frequency
    RGrep(tmp2, src1, "LABEL", "FN", false);
    RGrep(src1, src2, "NAME",  "CM");
    
    RFilterC(src1, src1_, std::set<Label> { "EXP_FREQ", "REF_DEPTH_NORMAL", "VAR_DEPTH_NORMAL", "OBS_FREQ_NORMAL", "REF_DEPTH_TUMOR", "VAR_DEPTH_TUMOR", "OBS_FREQ_TUMOR", "QUAL" }, true);
    RFilterC(src2, src2_, std::set<Label> { "EXP_FREQ", "REF_DEPTH_NORMAL", "VAR_DEPTH_NORMAL", "OBS_FREQ_NORMAL", "REF_DEPTH_TUMOR", "VAR_DEPTH_TUMOR", "OBS_FREQ_TUMOR", "QUAL" }, true);

    // Aggregated mean from all detected sequins
    RAggregateMean(src1_, tmpM, "EXP_FREQ", Imputation::Remove);
    
    // Aggregated standard deviation from only CM sequins
    RAggregateSD(src2_, tmpS, "EXP_FREQ", Imputation::Remove);

    auto read = [&](const Label &l, std::map<double, double> &m, std::map<double, double> &s)
    {
        // Average for each allele frequency group (all detected sequins)
        RFilterC(tmpM, tmp1, std::set<Label> { "EXP_FREQ", l}, true); m = RBinaryTSV(tmp1, "EXP_FREQ", l);

        // Standard deviation for each allele frequency group (only from CM)
        RFilterC(tmpS, tmp1, std::set<Label> { "EXP_FREQ", l}, true); s = RBinaryTSV(tmp1, "EXP_FREQ", l);
    };
    
    std::map<double, double> m1, m3, m4, m5, m6;
    std::map<double, double> s1, s3, s4, s5, s6;

    read("REF_DEPTH_NORMAL", m1, s1);
    read("REF_DEPTH_TUMOR",  m3, s3);
    read("VAR_DEPTH_TUMOR",  m4, s4);
    read("OBS_FREQ_TUMOR",   m5, s5);
    read("QUAL",             m6, s6);

    o.generate("somatic_variant_table.tsv");
    o.writer->open("somatic_variant_table.tsv");
    o.writer->write("EXP_FREQ\tOBS_FREQ\tTP\tFN\tREF_DEPTH_TUMOR\tVAR_DEPTH_TUMOR\tDEPTH_NORMAL\tQUAL_SCORE");

    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";
    
    // For each allele frequency group ...
    for (auto &i : m1)
    {
        auto ms = [&](const std::map<double, double> &m, const std::map<double, double> &s)
        {
            if (!m.count(i.first))
            {
                return std::string("NA");
            }
            else
            {
                return S4(m.at(i.first)) + " +- " + (s.count(i.first) ? S4(s.at(i.first)) : "NA");
            }
        };
        
        auto count = [&](const Label &l)
        {
            RGrep(src1_, tmp1, "EXP_FREQ", std::to_string(i.first), true, true);
            return RCount(tmp1, "EXP_FREQ", l);
        };

        o.writer->write((boost::format(f) % toString(i.first, 4)
                                          % ms(m5, s5)
                                          % count("TP")
                                          % count("FN")
                                          % ms(m3, s3)
                                          % ms(m4, s4)
                                          % ms(m1, s1)
                                          % ms(m6, s6)).str());
    }
    
    o.writer->close();
}

template <typename O> void writeROC(const FileName &file, const O &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotStrelkaROC()) % date()
                                                     % o.cmd
                                                     % o.work
                                                     % "somatic_sequin.tsv").str());
    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const FileName &es,
                         const FileName &ss,
                         const GSomatic::Stats &stats,
                         const GVariant::Options &o)
{
    const auto &r = Standard::instance().gen;
    const auto l1 = stats.oa.linear();
    
    Base   sl = 0;
    Counts sc = 0;
    
    for (auto &i : *(r.r2()))
    {
        for (auto &j : i.second.l2d)
        {
            if (isSoma(j.second.name) || isGerm(j.second.name) || isHP(j.second.name) || isMS(j.second.name))
            {
                sc++;
                sl += j.second.l.length();
            }
        }
    }
    
    auto S = [&](const GVariant::TableRow &x)
    {
        return !x.valid ? MISSING : (std::to_string(x.tp) + "/" + std::to_string(x.nr));
    };

    const auto f = "SEQUIN REPORT:                   %1%\n\n"
                   "REFERENCE FILES\n"
                   "Reference variants:              %2%\n"
                   "Reference regions:               %3%\n\n"
                   "USER-SUPPLIED FILES\n"
                   "User sample variants:            %4%\n"
                   "User sequin variants:            %5%\n"
                   "User regions:                    %6%\n\n"
                   "OUTPUT FILES\n"
                   "Somatic regions:                 %7%\n"
                   "Sample variants:                 %8%\n"
                   "True positive variants:          %9%\n"
                   "False positive variants:         %10%\n"
                   "False negative variants:         %11%\n"
                   "Variant summary:                 %12%\n"
                   "Summary report:                  %13%\n\n"
                   "GERMLINE VARIANTS\n"
                   "Sequin regions:                  %14%\n"
                   "Sequin size (bp):                %15%\n"
                   "Detected sample variants:        %16%\n"
                   "Detected sequin variants (TP):   %17%\n"
                   "Undetected sequin variants (FN): %18%\n"
                   "Erroneous sequin variants (FP):  %19%\n"
                   "False positives per KB:          %20%\n\n"
                   "SOMATIC VARIANT FREQUENCY:\n"
                   "Slope:                           %21%\n"
                   "R2:                              %22%\n"
                   "Somatic allele frequency table:  %23%";

        o.generate(file);
        o.writer->open(file);
        o.writer->write((boost::format(f) % date()
                                          % r.v1()->src
                                          % (o.combined ? (r.r1()->src + " and " + r.r3()->src) : r.r1()->src)
                                          % (es.empty() ? MISSING : es)
                                          % (ss.empty() ? MISSING : ss)
                                          % (o.uBED.empty() ? MISSING : o.uBED)
                                          % (o.work + "/genome_files/somatic_sequin.bed") // 7
                                          % (o.work + "/genome_files/somatic_sample.vcf") // 8
                                          % (o.work + "/genome_files/somatic_TP.vcf")     // 9
                                          % (o.work + "/genome_files/somatic_FP.vcf")     // 10
                                          % (o.work + "/genome_files/somatic_FN.vcf")     // 11
                                          % (o.work + "/somatic_sequin.tsv")              // 12
                                          % (o.work + "/somatic_summary.txt")             // 13
                                          % sc                 // 14
                                          % sl                 // 15
                                          % stats.hs.vs.size() // 16
                                          % stats.ds.oc.tp()   // 17
                                          % stats.ds.oc.fn()   // 18
                                          % stats.ds.oc.fp()   // 19
                                          % S2(1000.0 * ((float) stats.ds.oc.fp() / sl))
                                          % l1.m   // 21
                                          % l1.R2  // 22
                                          % (o.work + "/somatic_variant_table.tsv")).str());
    o.writer->close();
}

void GSomatic::report(const FileName &f1, const FileName &f2, const GVariant::Options &o)
{
    createD(o.work + "/genome_files");
    const auto &r = Standard::instance().gen;

    GSomatic::Stats stats;
    GVariant::analyze(stats, GSomatic(), f1, f2, o, Standard::instance().gen.v3());

    for (const auto &m : stats.ds.tps)
    {
        auto oAF = [&]()
        {
            /*
             * Since this tool works with both germline and somatic callers.
             * We need to guess what the observed allele frequency.
             */
            
            // Tumor allele fequency?
            const auto tAF = tumorAF(m.qry);
            
            return std::isnan(tAF) ? m.qry.allF : tAF;
        };

        const auto exp = r.af(m.var->name);
        const auto obs = oAF();
        
        A_ASSERT(!std::isnan(exp));
        
        // Eg: 2821292107
        const auto x = toString(m.var->key());
        
        // Add for all variants
        stats.oa.add(x, exp, obs);
        
        // Add for mutation type
        stats.m2a[m.qry.type()].add(x, exp, obs);
    }
    
    // Determine quantification limits
    stats.oa.limit = stats.oa.limitQuant();
    
    for (auto &i : stats.m2a)
    {
        i.second.limit = i.second.limitQuant();
    }
    
    /*
     * Performance by allele frequency
     */
    
    for (const auto &i : r.v1()->data.vars())
    {
        if (isSoma(i.name))
        {
            const auto af = r.af(i.name);
            stats.af[af].nr()++;
            
            if (stats.ds.findTP(i.name))
            {
                stats.af[af].tp()++;
            }
            else
            {
                stats.af[af].fn()++;
            }
        }
    }

    // Generating somatic_sequin.tsv
    writeSomatic("somatic_sequin.tsv", stats, o);

    // Generating somatic_FN.vcf
    writeFN("somatic_FN.vcf", stats.ds.fns, o);

    // Generating somatic_features.bed
    copy(Standard::instance().gen.a1()->src, o.work + "/genome_files/somatic_features.bed");

    // Generating somatic_variant_table.tsv
    writeTable("somatic_variant_table.tsv", stats, o);

    // Generating somatic_summary.txt
    writeSummary("somatic_summary.txt", f1, f2, stats, o);

    if (o.report)
    {
        // Generating somatic_ROC.R
        writeROC("report_files/somatic_ROC.R", o);
        
        // Generating somatic_qualFilter.R
        writeQualFilter(stats, o);
        
        // Generating somatic_ladder.R
        writeAllele(o);
        
        Report::somatic(o);
    }
}
