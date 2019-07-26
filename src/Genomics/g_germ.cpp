#include <cmath>
#include "tools/report.hpp"
#include "data/resources.hpp"
#include "Genomics/g_germ.hpp"
#include "writers/vcf_writer.hpp"

using namespace Anaquin;

typedef GVariant::Options Options;

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
        
        const auto name = (i.var && i.alt && i.ref ? i.var->name : MISSING);
        o.writer->write((boost::format(format) % i.rID
                                               % i.qry.cID
                                               % i.qry.l.start
                                               % label
                                               % gt2str(i.qry.gt)
                                               % var2str(i.qry.type())
                                               % (name != MISSING ? std::to_string(r.af(name)) : MISSING)
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

template <typename T, typename O> void writeSV(const std::string &format,
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
    o.writer->write((boost::format(format) % "NAME"
                                           % "CHROM"
                                           % "POSITION"
                                           % "LABEL"
                                           % "GENOTYPE"
                                           % "TYPE"
                                           % "EXP_FREQ"
                                           % "REF_DEPTH"
                                           % "VAR_DEPTH"
                                           % "OBS_FREQ"
                                           % "QUAL"
                                           % "SIZE"
                                           % r.a1()->strForKeys()).str());

    writeGQ(format, stats, o); // Write true positives and false negatives
    writeGD(format, stats, o); // Write false positives
    writeSV(format, stats, o); // Write sample variants
    
    o.writer->close();
}

static void writeTable(const FileName &file, const GGerm::Stats &stats, const GVariant::Options &o)
{
    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write("VARIANTS\tTOTAL\tTP\tFP\tFN\tSN\tPC\tFP_KB\tDEPTH\tSAMPLE");

    auto write = [&](const Label &name, const GVariant::TableRow &x)
    {
        o.writer->write((boost::format(f) % name
                                          % (x.valid ? replaceNA(x.nr)     : MISSING)
                                          % (x.valid ? replaceNA(x.tp)     : MISSING)
                                          % (x.valid ? replaceNA(x.fp)     : MISSING)
                                          % (x.valid ? replaceNA(x.fn)     : MISSING)
                                          % (x.valid ? replaceNA(x.sn)     : MISSING)
                                          % (x.valid ? replaceNA(x.pc)     : MISSING)
                                          % (x.valid ? replaceNA(x.fpKB)   : MISSING)
                                          % (x.valid ? replaceNA(x.depth)  : MISSING)
                                          % (x.valid ? replaceNA(x.sample) : MISSING)).str());
    };

    write("Total_variants",              GVariant::getTotal(o));
    write("SNV_type",                    GVariant::getSNV(o));
    write("Indel_type",                  GVariant::getIndel(o));
    write("Homozygous_Genotype",         GVariant::getHom(o));
    write("Heterozygous_Genotype",       GVariant::getHet(o));
    write("CodingRegion_GeneContext",    GVariant::getCode(o));
    write("NoncodingRegion_GeneContext", GVariant::getNCode(o));
    write("ATrich_GCcontent",            GVariant::getAT(o));
    write("GCrich_GCcontent",            GVariant::getGC(o));
    write("SINE_MobileElement",          GVariant::getSine(o));
    write("LINE_MobileElement",          GVariant::getLine(o));
    write("LTR_MobileElement",           GVariant::getLTR(o));
    write("DNA_repeat_MobileElement",    GVariant::getDNA(o));
    write("Mono_SimpleRepeat",           GVariant::getDi(o));
    write("Di_SimpleRepeat",             GVariant::getTri(o));
    write("Tri_SimpleRepeat",            GVariant::getMono(o));
    write("Quad_SimpleRepeat",           GVariant::getQuad(o));

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
            {
                sc++;
                sl += j.second.l.length();
            }
        }
    }

    const auto f = "SEQUIN REPORT:                   %1%\n\n"
                   "REFERENCE FILES\n"
                   "Reference variants:              %2%\n"
                   "Reference regions:               %3%\n\n"
                   "USER-SUPPLIED FILES\n"
                   "User sample variants:            %4%\n"
                   "User sequin variants:            %5%\n"
                   "User regions:                    %6%\n\n"
                   "OUTPUT FILES\n"
                   "Germline regions:                %7%\n"
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
                   "DIAGNOSTIC PERFORMANCE\n"
                   "Total variants:                  %21%\n"
                   "SNV (Mutation):                  %22%\n"
                   "Indel (Mutation):                %23%\n"
                   "Homozygous (Genotype):           %24%\n"
                   "Heterozygous (Genotype):         %25%\n"
                   "Coding region (GeneContext):     %26%\n"
                   "Noncoding region (GeneContext):  %27%\n"
                   "ATrich (GCcontent):              %28%\n"
                   "GCrich (GCcontent):              %29%\n"
                   "SINE (MobileElement):            %30%\n"
                   "DNA repeat (MobileElement):      %31%\n"
                   "LINE (MobileElement):            %32%\n"
                   "LTR (MobileElement):             %33%\n"
                   "Mono (Simple repeat):            %34%\n"
                   "Di (Simple repeat):              %35%\n"
                   "Tri (Simple repeat):             %36%\n"
                   "Quad (Simple repeat):            %37%\n"
                   "Germline variant table:          %38%";
    
    const auto rr = o.combined ? (r.r1()->src + " and " + r.r3()->src) : r.r1()->src;
    
    auto S = [&](const GVariant::TableRow &x)
    {
        return !x.valid ? MISSING : (std::to_string(x.tp) + "/" + std::to_string(x.nr));
    };

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % date()
                                      % r.v1()->src
                                      % rr
                                      % (es.empty() ? MISSING : es)
                                      % (ss.empty() ? MISSING : ss)
                                      % (o.uBED.empty() ? MISSING : o.uBED)
                                      % (o.work + "/genome_files/germline_sequin.bed") // 7
                                      % (o.work + "/genome_files/germline_sample.vcf") // 8
                                      % (o.work + "/genome_files/germline_TP.vcf")     // 9
                                      % (o.work + "/genome_files/germline_FP.vcf")     // 10
                                      % (o.work + "/genome_files/germline_FN.vcf")     // 11
                                      % (o.work + "/germline_sequin.tsv")              // 12
                                      % (o.work + "/germline_summary.txt")             // 13
                                      % sc                 // 14
                                      % sl                 // 15
                                      % stats.hs.vs.size() // 16
                                      % stats.ds.oc.tp()   // 17
                                      % stats.ds.oc.fn()   // 18
                                      % stats.ds.oc.fp()   // 19
                                      % S2(1000.0 * ((float) stats.ds.oc.fp() / sl))
                                      % stats.ds.oc.nr()         // 21
                                      % S(GVariant::getSNV(o))   // 22
                                      % S(GVariant::getIndel(o)) // 23
                                      % S(GVariant::getHom(o))   // 24
                                      % S(GVariant::getHet(o))   // 25
                                      % S(GVariant::getCode(o))  // 26
                                      % S(GVariant::getNCode(o)) // 27
                                      % S(GVariant::getAT(o))    // 28
                                      % S(GVariant::getGC(o))    // 29
                                      % S(GVariant::getSine(o))  // 30
                                      % S(GVariant::getDNA(o))   // 31
                                      % S(GVariant::getLine(o))  // 32
                                      % S(GVariant::getLTR(o))   // 33
                                      % S(GVariant::getMono(o))  // 34
                                      % S(GVariant::getDi(o))    // 35
                                      % S(GVariant::getTri(o))   // 36
                                      % S(GVariant::getQuad(o))  // 37
                                      % (o.work + "/germline_variant_table.tsv")
                     ).str());
    o.writer->close();
}

void GGerm::report(const FileName &f1, const FileName &f2, const Options &o)
{
    createD(o.work + "/genome_files");

    GVariant::Stats stats;
    GVariant::analyze(stats, GGerm(), f1, f2, o, Standard::instance().gen.v2());
    
    // Generating germline_sequin.tsv
    writeQuins("germline_sequin.tsv", stats, o);
    
    // Generating germline_FN.vcf
    writeFN("germline_FN.vcf", stats.ds.fns, o);
    
    // Generating germline_features.bed
    copy(Standard::instance().gen.a1()->src, o.work + "/genome_files/germline_features.bed");

    // Generating germline_variant_table.tsv
    writeTable("germline_variant_table.tsv", stats, o);
    
    // Generating germline_summary.txt
    writeSummary("germline_summary.txt", f1, f2, stats, o);

    if (o.report)
    {
        Report::germline(o);
    }
}
