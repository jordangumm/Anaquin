#include <cmath>
#include "tools/report.hpp"
#include "Genomics/g_somatic.hpp"

using namespace Anaquin;

template <typename T, typename O> void writeQualFilter(const T &, const O &o)
{
    o.generate("somatic_qualFilter.R");
    o.writer->open("genome_files/" + o.name + "_qualFilter.R");
    o.writer->write((boost::format(PlotQualFilter()) % date()
                                                     % o.cmd
                                                     % o.work
                                                     % (o.name + "_sequin.tsv")).str());
    o.writer->close();
}

template <typename O> void writeAllele(const O &o)
{
    o.generate("somatic_ladder.R");
    o.writer->open("genome_files/" + o.name + "_ladder.R");
    o.writer->write((boost::format(PlotAllele()) % date()
                                                 % o.cmd
                                                 % o.work
                                                 % (o.name + "_sequin.tsv")).str());
    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const FileName &es,
                         const FileName &ss,
                         const GSomatic::Stats &stats,
                         const GVariant::Options &o)
{
    std::stringstream tab;

    for (auto iter = stats.af.rbegin(); iter != stats.af.rend(); iter++)
    {
        const auto format = "       %1$.5f%2%%3%/%4%\n";
        
        const auto tp = iter->second.tp();
        const auto fn = iter->second.fn();
        
        tab << (boost::format(format) % iter->first
                                      % std::string(16, ' ')
                                      % tp
                                      % (tp + fn)).str();
    }
    
    const auto summary = "-------SUMMARY STATISTICS\n\n"
                         "-------REFERENCE FILES\n\n"
                         "       Reference variants: %1%\n"
                         "       %2%\n\n"
                         "-------USER-SUPPLIED FILES\n\n"
                         "       Input file (first):  %3%\n"
                         "       Input file (second): %4%\n"
                         "       User regions:        %5%\n\n"
                         "-------OUTPUT FILES\n\n"
                         "       Somatic regions:         %6%\n"
                         "       True positive variants:  %7%\n"
                         "       False positive variants: %8%\n"
                         "       False negative variants: %9%\n"
                         "       Sample variants:         %10%\n"
                         "       Variant summary:         %11%\n"
                         "       Summary report:          %12%\n\n"
                         "-------SOMATIC VARIANTS\n\n"
                         "       Sequin regions:                  %13%\n"
                         "       Sequin size (bp):                %14%\n"
                         "       Detected sample variants:        %15%\n"
                         "       Detected sequin variants (TP):   %16%\n"
                         "       Undetected sequin variants (FN): %17%\n"
                         "       Erroneous sequin variants (FP):  %18%\n"
                         "       False positives per kb:          %19$.5f\n\n"
                         "-------LINEAR REGRESSION\n\n"
                         "       Slope: %20%\n"
                         "       R2:    %21%\n\n"
                         "-------DIAGNOSTIC PERFORMANCE\n\n"
                         "       EXPECTED ALLELE FQ     DETECTED VARIANTS\n"
                         "%22%";

        const auto &r = Standard::instance().gen;
        const auto lm = stats.oa.linear();
    
        Base   sl = 0;
        Counts sc = 0;
    
        for (auto &i : *(r.r2()))
        {
            for (auto &j : i.second.l2d)
            {
                if (isSoma(j.second.name) || isGerm(j.second.name) || isHP(j.second.name) || isMS(j.second.name))
                { sc++; sl += j.second.l.length(); }
            }
        }

        auto rr = o.combined ? ("Reference regions (human): " + r.r1()->src + "\n" + "       Reference regions (decoy): " + r.r3()->src) :
                                "Reference regions: " + r.r1()->src;

        o.generate(file);
        o.writer->open(file);
        o.writer->write((boost::format(summary) % r.v1()->src
                                                % rr
                                                % (es.empty() ? "-" : es)
                                                % (ss.empty() ? "-" : ss)
                                                % (o.uBED.empty() ? "-" : o.uBED)
                                                % (o.work + "/genome_files/somatic_sequin.bed")
                                                % (o.work + "/genome_files/somatic_TP.vcf")
                                                % (o.work + "/genome_files/somatic_FP.vcf")
                                                % (o.work + "/genome_files/somatic_FN.vcf")
                                                % (o.work + "/genome_files/somatic_sample.vcf")
                                                % (o.work + "/somatic_sequin.tsv")
                                                % (o.work + "/somatic_summary.stats")
                                                % sc // 13
                                                % sl // 14
                                                % stats.hs.vs.size()
                                                % stats.ds.oc.tp() // 16
                                                % stats.ds.oc.fn() // 17
                                                % stats.ds.oc.fp() // 18
                                                % (1000.0 * ((float) stats.ds.oc.fp() / sl))
                                                % lm.m
                                                % lm.R2
                                                % tab.str()).str());
    o.writer->close();
}

void GSomatic::report(const FileName &f1, const FileName &f2, const GVariant::Options &o)
{
    const auto &r = Standard::instance().gen;
    
    createD(o.work + "/genome_files");
    createD(o.work + "/report_files");

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

    // Generating somatic_ladder.R
    writeAllele(o);

    // Generating somatic_FN.vcf
    writeFN("somatic_FN.vcf", stats.ds.fns, o);

    // Generating somatic_ROC.R
    writeROC("somatic_ROC.R", o);

    // Generating somatic_qualFilter.R
    writeQualFilter(stats, o);

    // Generating somatic_features.bed
    copy(Standard::instance().gen.a1()->src, o.work + "/genome_files/somatic_features.bed");

    // Generating somatic_summary.stats
    writeSummary("somatic_summary.stats", f1, f2, stats, o);
}
