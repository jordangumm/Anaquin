#include <cmath>
#include "tools/report.hpp"
#include "data/resources.hpp"
#include "Genomics/g_germ.hpp"
#include "writers/vcf_writer.hpp"

using namespace Anaquin;

template <typename T, typename O> void writeFP(const std::vector<T> &x,
                                               const std::string &label,
                                               const std::string &format,
                                               const O &o)
{
    const auto &r = Standard::instance().gen;
    
    for (const auto &i : x)
    {
        if (i.rID.empty())
        {
            continue;
        }
        
        const auto sID = (i.var && i.alt && i.ref ? i.var->name : MISSING);
        o.writer->write((boost::format(format) % i.rID
                                               % i.qry.cID
                                               % i.qry.l.start
                                               % label
                                               % gt2str(i.qry.gt)
                                               % var2str(i.qry.type())
                                               % (sID != "-" ? std::to_string(r.af(sID)) : MISSING)
                                               % i.qry.dp(0)
                                               % i.qry.dp(1)
                                               % i.qry.obsAF()
                                               % toString(i.qry.qual[0])
                                               % MISSING
                                               % r.a1()->strForLocus(i.qry.cID, i.qry.l)).str());
    }
}

template <typename T, typename O> void writeGD(const std::string &format,
                                               const T &x,
                                               const O &o)
{
    writeFP(x.ds.fps, "FP", format, o);
}

template <typename T, typename O> void writeGS(const std::string &format,
                                               const T &x,
                                               const O &o)
{
    const auto &r = Standard::instance().gen;
    
    for (const auto &i : x.hs.vs)
    {
        o.writer->write((boost::format(format) % i.name
                                               % i.cID
                                               % i.l.start
                                               % "SV"
                                               % gt2str(i.gt)
                                               % var2str(i.type())
                                               % MISSING
                                               % i.dp(0)
                                               % i.dp(1)
                                               % i.obsAF()
                                               % toString(i.qual[0])
                                               % MISSING
                                               % r.a1()->strForLocus(i.cID, i.l)).str());
    }
}

template <typename T, typename O> void writeGQ(const std::string &format,
                                               const T &stats,
                                               const O &o)
{
    const auto &r = Standard::instance().gen;
    const auto r4 = r.r4()->ginters();
    const auto r5 = r.r5()->ginters();

    for (const auto &i : r.v1()->data.vars())
    {
        if (isGerm(i.name))
        {
            // Can we find this sequin?
            const auto isTP = stats.ds.findTP(i.name);
            
            // Locus for the sequin
            const auto l5 = r5.at(i.cID).overlap(i.l)->l();
            
            // This shouldn't fail...
            const auto &v = r.v1()->data.findByKey(i.key());

            if (isTP)
            {
                // Called variant (if found)
                const auto &q = isTP->qry;

                o.writer->write((boost::format(format) % i.name
                                                       % i.cID
                                                       % i.l.start
                                                       % "TP"
                                                       % gt2str(v->gt)
                                                       % var2str(i.type())
                                                       % r.af(i.name)
                                                       % q.dp(0)
                                                       % q.dp(1)
                                                       % q.obsAF()
                                                       % toString(q.qual[0])
                                                       % r4.at(i.cID).length(l5)
                                                       % r.a1()->strForLocus(i.cID, i.l)).str());
            }

            // Failed to detect the variant
            else
            {
                o.writer->write((boost::format(format) % i.name
                                                       % i.cID
                                                       % i.l.start
                                                       % "FN"
                                                       % gt2str(v->gt)
                                                       % var2str(i.type())
                                                       % r.af(i.name)
                                                       % MISSING
                                                       % MISSING
                                                       % MISSING
                                                       % MISSING
                                                       % r4.at(i.cID).length(l5)
                                                       % r.a1()->strForLocus(i.cID, i.l)).str());
            }
        }
    }
}

template <typename T, typename O> void writeQuins(const FileName &file, const T &stats, const O &o)
{
    const auto &r = Standard::instance().gen;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%%13%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Position"
                                           % "Label"
                                           % "Genotype"
                                           % "Type"
                                           % "ExpFreq"
                                           % "Ref.Depth"
                                           % "Var.Depth"
                                           % "Obs.Freq"
                                           % "Qual"
                                           % "Size"
                                           % r.a1()->strForKeys()).str());

    writeGQ(format, stats, o);
    writeGD(format, stats, o);
    writeGS(format, stats, o);
    
    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const FileName &es,
                         const FileName &ss,
                         const GGerm::Stats &stats,
                         const GVariant::Options &o)
{
    const auto &r = Standard::instance().gen;
    
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

    const auto summary = "-------SUMMARY STATISTICS\n\n"
                         "-------REFERENCE FILES\n\n"
                         "       Reference variants: %1%\n"
                         "       %2%\n\n"
                         "-------USER-SUPPLIED FILES\n\n"
                         "       Input file (first):      %3%\n"
                         "       Input file (second):     %4%\n"
                         "       User regions:            %5%\n\n"
                         "-------OUTPUT FILES\n\n"
                         "       Germline regions:        %6%\n"
                         "       True positive variants:  %7%\n"
                         "       False positive variants: %8%\n"
                         "       False negative variants: %9%\n"
                         "       Sample variants:         %10%\n"
                         "       Variant summary:         %11%\n"
                         "       Summary report:          %12%\n\n"
                         "-------GERMLINE VARIANTS\n\n"
                         "       Sequin regions:                  %13%\n"
                         "       Sequin size (bp):                %14%\n"
                         "       Detected sample variants:        %15%\n"
                         "       Detected sequin variants (TP):   %16%\n"
                         "       Undetected sequin variants (FN): %17%\n"
                         "       Erroneous sequin variants (FP):  %18%\n"
                         "       False positives per kb:          %19$.5f\n";

    auto rr = o.combined ? ("Reference regions (human): " + r.r1()->src + "\n" + "       Reference regions (decoy): " + r.r3()->src) :
                            "Reference regions: " + r.r1()->src;
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(summary) % r.v1()->src
                                            % rr
                                            % (es.empty() ? "-" : es)
                                            % (ss.empty() ? "-" : ss)
                                            % (o.uBED.empty() ? "-" : o.uBED)
                                            % (o.work + "/genome_files/germline_sequin.bed")
                                            % (o.work + "/genome_files/germline_TP.vcf")
                                            % (o.work + "/genome_files/germline_FP.vcf")
                                            % (o.work + "/genome_files/germline_FN.vcf")
                                            % (o.work + "/genome_files/germline_sample.vcf")
                                            % (o.work + "/germline_sequin.tsv")
                                            % (o.work + "/germline_summary.stats")
                                            % sc // 13
                                            % sl // 14
                                            % stats.hs.vs.size()
                                            % stats.ds.oc.tp() // 16
                                            % stats.ds.oc.fn() // 17
                                            % stats.ds.oc.fp() // 18
                                            % (1000.0 * ((float) stats.ds.oc.fp() / sl))).str());
    o.writer->close();
}

void GGerm::report(const FileName &f1, const FileName &f2, const Options &o)
{
    createD(o.work + "/genome_files");
    createD(o.work + "/report_files");

    GVariant::Stats stats;
    GVariant::analyze(stats, GGerm(), f1, f2, o, Standard::instance().gen.v2());
    
    // Generating germline_sequin.tsv
    writeQuins("germline_sequin.tsv", stats, o);
    
    // Generating germline_FN.vcf
    writeFN("germline_FN.vcf", stats.ds.fns, o);
    
    // Generating germline_features.bed
    copy(Standard::instance().gen.a1()->src, o.work + "/genome_files/germline_features.bed");

    // Generating germline_summary.stats
    writeSummary("germline_summary.stats", f1, f2, stats, o);
}
