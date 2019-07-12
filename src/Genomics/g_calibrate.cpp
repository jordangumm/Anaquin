#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include "tools/report.hpp"
#include "tools/random.hpp"
#include "Genomics/Genomics.hpp"
#include "writers/bam_writer.hpp"
#include "Genomics/g_calibrate.hpp"

using namespace Anaquin;

// Defined in ProcessReads
extern Anaquin::AlignedReads __CombinedBAM1__;

// Defined in ProcessReads
extern bam_hdr_t *__CombinedBAMHeader__;

// Defined in ProcessReads
extern std::vector<const Anaquin::ReadName *> __HackBAMDecoy2__;

typedef ParserBAMBED::Stats Stats;
typedef GCalibrate::Method Method;
typedef GCalibrate::Options Options;
typedef ParserBAMBED::Response Response;

static ParserBAMBED::Stats parseSamp(const FileName &file, const Chr2DInters &r, GCalibrate::CalibrateStats &stats, const GCalibrate::Options &o)
{
    return ParserBAMBED::parse(file, r, [&](const ParserBAM::Data &x, const ParserBAM::Info &info, const DInter *)
    {
        if (info.p && !(info.p % 1000000)) { o.wait(toString(info.p)); }
        
        if (x.mapped)
        {
            stats.nh++;
        }
        
        return ParserBAMBED::Response::OK;
    });
}

static ParserBAMBED::Stats parseSeq(const FileName &file, const Chr2DInters &r, GCalibrate::CalibrateStats &stats, const GCalibrate::Options &o)
{
    return ParserBAMBED::parse(file, r, [&](ParserBAM::Data &x, const ParserBAM::Info &info, const DInter *)
    {
        if (info.p && !(info.p % 1000000)) { o.wait(toString(info.p)); }
        
        if (x.mapped)
        {
            stats.nd++;
        }
        
        x.lName();
        return stats.trimmed.count(x.name) ? Response::SKIP_EVERYTHING : Response::OK;
    });
}

Stats GCalibrate::calibrate(const FileName &src,
                            const FileName &dst,
                            const NormFactors &norms,
                            const Chr2DInters &r1,
                            const Chr2DInters &r2,
                            const std::set<ReadName> &trimmed,
                            CalibrateStats &,
                            const Options &o)
{
    typedef std::map<ChrID, std::map<Locus, std::shared_ptr<RandomSelection>>> Selection;
    
    Selection select;
    
    for (const auto &i : norms)
    {
        for (const auto &j : i.second)
        {
            A_ASSERT(j.second >= 0 && j.second <= 1.0 && !std::isnan(j.second));
            
            // Create independent random generator for each calibration region
            select[i.first][j.first] = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - j.second));
        }
    }

    A_ASSERT(select.size() == norms.size());

    BAMWriter w;
    w.open(dst);
    
    auto after = ParserBAMBED::parse(src, r1, [&](ParserBAM::Data &x, const ParserBAM::Info &i, const DInter *)
    {
        if (i.p && !(i.p % 1000000)) { o.logInfo(toString(i.p)); }

        if (i.p == 0)
        {
            // Write headers in case there's no sample read
            w.writeH(x);
        }
        
        x.lName();
        
        if (trimmed.count(x.name))
        {
            return Response::SKIP_EVERYTHING;
        }
        
        auto kept = !x.mapped || x.isDuplicate;
        
        if (!kept)
        {
            DInter *inter = nullptr;
            
            /*
             * Should that be contains or overlap? We prefer overlaps because any read that is overlapped
             * into the regions still give valuable information and sequencing depth.
             */

            if (r1.count(x.cID) && (inter = r1.at(x.cID).overlap(x.l)))
            {
                // Calibrate both paired-ends because they have the same name
                if (select.at(x.cID).at(inter->l())->select(x.name))
                {
                    kept = true;
                }
            }
            else
            {
                // Never throw away reads outside the regions
                kept = true;
            }
        }
        
        if (kept)
        {
            w.write(x);
            
            if (r2.count(x.cID) && r2.at(x.cID).overlap(x.l) && !x.isDuplicate)
            {
                r2.at(x.cID).overlap(x.l)->map(x.l);
            }
            
            return Response::OK;
        }
        
        return Response::SKIP_EVERYTHING;
    });

    w.close();
    return after;
}

template <typename Stats> Coverage stats2cov(const GCalibrate::Method meth, const Stats &stats)
{
    switch (meth)
    {
        case Method::Mean:
        case Method::Percent: { return stats.mean; }
        case Method::Median:  { return stats.p50;  }
    }
}

static bool shouldTrim(const Locus &l, const DInter *m, const GCalibrate::Options &o)
{
    const auto lTrim = std::abs(l.start - m->l().start) <= o.trim;
    const auto rTrim = std::abs(l.end - m->l().end) <= o.trim;
    return lTrim || rTrim;
}

static void trimming(const FileName &file, const Chr2DInters &r1, GCalibrate::CalibrateStats &stats, const GCalibrate::Options &o)
{
    if (o.trim)
    {
        std::vector<DInter *> multi;
        o.logInfo("Checking trimming");
        
        ParserBAMBED::parse(file, r1, [&](ParserBAM::Data &x, const ParserBAM::Info &info, const DInter *)
        {
            if (info.p && !(info.p % 1000000)) { o.wait(std::to_string(info.p)); }
            
            multi.clear();
            const auto m = x.mapped && r1.count(x.cID) ? r1.at(x.cID).overlap(x.l, &multi) : nullptr;
            
            if (m)
            {
                std::sort(multi.begin(), multi.end(), [&](const DInter * x, const DInter * y)
                {
                    return x->l().length() < y->l().length();
                });
                
                for (const auto &i : multi)
                {
                    if (shouldTrim(x.l, i, o))
                    {
                        x.lName();
                        stats.trimmed.insert(x.name);
                        break;
                    }
                }
            }
            
            return ParserBAMBED::Response::OK;
        });
    }
}

static void calcNorm(const Chr2DInters &r1, const Chr2DInters &r2, GCalibrate::CalibrateStats &stats, const GCalibrate::Options &o)
{
    // Cache index by sequin ID, more convenient to work with
    std::map<SequinID, const DInter *> cES, cSS;
    
    for (auto &i : stats.es.inters)
    {
        for (auto &j : i.second.data()) { cES[j.first] = &j.second; }
    }

    for (auto &i : stats.ss.inters)
    {
        for (auto &j : i.second.data()) { cSS[j.first] = &j.second; }
    }

    // For each chromosome...
    for (auto &i : r1)
    {
        const auto &cID = i.first;
        
        // For each region...
        for (auto &j1 : i.second.data())
        {
            const auto &l1 = j1.second.l();
            
            /*
             * We'll need to sample all reads in r1 but we might have multiple calibration regions in r2.
             * Let's combine them.
             */
            
            std::vector<DInter *> o2;
            
            if (r2.count(cID))
            {
                r2.at(cID).overlap(l1, &o2);
            }
            
            DInter::Stats g1, s1;
            
            if (o2.empty())
            {
                g1.sums = g1.length = g1.mean = g1.p50 = 0;
                s1.sums = s1.length = s1.mean = s1.p50 = 0;
            }
            else
            {
                std::vector<DInter::Stats> g2, s2;
                
                for (const auto &j2 : o2)
                {
                    const auto &l = j2->l();
                    
                    switch (o.cMode)
                    {
                        case GCalibrate::CalibrateMode::TwoBAM:
                        {
                            g2.push_back(stats.es.inters.at(cID).find(l.key())->stats());
                            s2.push_back(stats.ss.inters.at(cID).find(l.key())->stats());
                            break;
                        }

                        case GCalibrate::CalibrateMode::Combined:
                        {
                            g2.push_back(cES[j2->name()]->stats());
                            s2.push_back(cSS[j2->name()]->stats());
                            break;
                        }
                    }
                }

                assert(!g2.empty() && !s2.empty());
                
                for (const auto &i : g2)
                {
                    g1.sums += i.sums; g1.aligns += i.aligns; g1.length += i.length; g1.nonZeros += i.nonZeros; g1.zeros += i.zeros;
                    for (const auto &j : i.raws) { g1.raws.push_back(j); }
                }
                
                for (const auto &i : s2)
                {
                    s1.sums += i.sums; s1.aligns += i.aligns; s1.length += i.length; s1.nonZeros += i.nonZeros; s1.zeros += i.zeros;
                    for (const auto &j : i.raws) { s1.raws.push_back(j); }
                }
                
                std::sort(g1.raws.begin(), g1.raws.end());
                std::sort(s1.raws.begin(), s1.raws.end());
                
                g1.p50  = SS::quant(g1.raws, 0.50); // Combined statistics
                s1.p50  = SS::quant(s1.raws, 0.50); // Combined statistics
                g1.mean = g1.sums / g1.length;      // Combined statistics
                s1.mean = s1.sums / s1.length;      // Combined statistics
            }
            
            /*
             * Now we have the data, we'll need to compare coverage and determine the fraction that
             * the synthetic alignments needs to be sampled.
             */
            
            const auto seqsC = stats2cov(o.meth, s1);
            const auto endoC = o.meth == Method::Percent ? seqsC : stats2cov(o.meth, g1);
            
            Proportion norm;
            
            switch (o.meth)
            {
                case Method::Mean:
                case Method::Median:  { norm = seqsC ? std::min(endoC / seqsC, 1.0) : NAN; break; }
                case Method::Percent: { norm = o.p; break; }
            }
            
            stats.allBeforeEndoC.push_back(endoC);
            stats.allBeforeSeqsC.push_back(seqsC);
            
            if (std::isnan(norm))
            {
                o.logWarn((boost::format("Normalization is zero for %1%:%2%") % i.first % l1.start).str());
                
                // We can't use NAN...
                norm = 0.0;
            }
            else if (norm == 1.0)
            {
                o.logWarn((boost::format("Normalization is 1 for %1%:%2%-%3% (%4%)") % i.first
                                                                                     % l1.start
                                                                                     % l1.end
                                                                                     % j1.first).str());
            }
            
            stats.c2v[cID][l1].nEndo   = g1.aligns;
            stats.c2v[cID][l1].nBefore = s1.aligns;
            
            stats.c2v[cID][l1].rID    = j1.first;
            stats.c2v[cID][l1].endo   = endoC;
            stats.c2v[cID][l1].before = seqsC;
            stats.c2v[cID][l1].norm   = stats.norms[i.first][l1] = norm;
            
            stats.allNorms.push_back(norm);
        }
    }
}

GCalibrate::CalibrateStats GCalibrate::combined(const FileName &file,
                                                const FileName &tmp,
                                                const Chr2DInters &r2,
                                                const Chr2DInters &r3,
                                                const Chr2DInters &r4,
                                                SStats &ss,
                                                const Options &o)
{
    GCalibrate::CalibrateStats stats;

    stats.es.inters = r3;
    stats.ss.inters = r4; // Calibration is only done in trimmed regions
    
    for (auto &i : stats.es.inters) { i.second.build(); }
    for (auto &i : stats.ss.inters) { i.second.build(); }

    BAMWriter w1;
    w1.open(tmp);

    BAMWriter w2;
    w2.open(o.work + "/calibrate_sample.bam");
    
    __CombinedBAM1__.clear(); __HackBAMDecoy2__.clear();

    ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        if (i.p && !(i.p % 1000000)) { o.wait(toString(i.p)); }

        if (i.p == 0)
        {
            // Required when constructing the header later...
            __CombinedBAMHeader__ = (bam_hdr_t *) x.copyH();
            
            if (!o.onlySeqLad)
            {
                // Write headers in case there's no sample read
                w2.writeH(x);
            }
        }

        // Always print out unmapped reads
        if (!x.mapped) { w2.write(x); return; }

        // Overlapped with mirror human regions?
        DInter *hr = stats.es.inters.count(x.cID) ? stats.es.inters.at(x.cID).overlap(x.l) : nullptr;
        
        // Overlapped with non-trimmed decoy regions?
        const DInter *s2 = !hr && r2.count(x.cID) ? r2.at(x.cID).overlap(x.l) : nullptr;
        
        // Overlapped with trimmed decoy regions?
        DInter *s4 = !hr && stats.ss.inters.count(x.cID) ? stats.ss.inters.at(x.cID).overlap(x.l) : nullptr;
        
        const auto decoy = isDecoy(x.cID);
        
        if (!decoy) { stats.nh++; }
        else        { stats.nd++; }
        
        /*
         * We don't assume sorted BAM, and there is no way to determine if both paired-end reads
         * are present. We will only count on the first pair.
         */
        
        if (!decoy && x.isFirstPair)
        {
            ss.hackHumanReads++;
        }

        if (hr)
        {
            if (x.mapped) { stats.es.nMap++; }
            else          { stats.es.nNA++;  }
            
            hr->map(x.l);
            assert(!s2 && !s4);
        }
        
        if (s2)
        {
            if (x.mapped) { stats.ss.nMap++; }
            else          { stats.ss.nNA++;  }

            if (!shouldTrim(x.l, s2, o))
            {
                if (s4) { s4->map(x.l); }
                w1.write(x);
            }
        }
        else if (isBegin(x.cID, "chrQL"))
        {
            w1.write(x); // Not calibrated but not sample reads
        }
        else
        {
            // Always print out anything that is not in the calibrated regions
            w2.write(x);
            
            // Decoy chromosome (e.g. chrQS) but outside calibrated regions?
            if (decoy)
            {
                stats.nSeqOut++;
            }
        }
    });

    w1.close();
    w2.close();
    
    return stats;
}

GCalibrate::CalibrateStats GCalibrate::twoBAM(const FileName &endo,
                                              const FileName &seqs,
                                              const Chr2DInters &r1,
                                              const Chr2DInters &r2,
                                              const Options &o)
{
    GCalibrate::CalibrateStats stats;
    
    // Check sample alignments before calibration
    stats.es = parseSamp(endo, r2, stats, o);
    
    // Calculate trimming on sequin alignments
    trimming(seqs, r1, stats, o);
    
    // Check sequin alignments before calibration
    stats.ss = parseSeq(seqs, r2, stats, o);
    
    return stats;
}

GCalibrate::Stats GCalibrate::analyze(const FileName &f1, const FileName &f2, const Options &o)
{
    assert(o.meth != Method::Percent || !std::isnan(o.p));
    const auto &r = Standard::instance().gen;
    
    o.analyze(f1);
    if (!f2.empty()) { o.analyze(f2); }
    
    o.logInfo("Trim: " + std::string(o.trim ? "true" : "false"));
    o.logInfo("Edge: " + toString(o.edge, 0));
    
    GCalibrate::Stats stats;
    
    /*
     * Combined: Human regions no trimming
     *           Decoy regions no trimming
     *           Human regions trimmed
     *           Decoy regions trimmed
     *
     * Both BAM: Human + decoy regions no trimming
     *           Human + decoy regions trimmed
     */
    
    const auto r1 = r.r1()->inters();
    const auto r2 = r.r2()->inters();
    const auto r3 = r.r3() ? r.r3()->inters() : r1;
    const auto r4 = r.r4() ? r.r4()->inters() : r2;

    // BAM file for combined mode
    const auto tmp = o.cMode == CalibrateMode::Combined ? tmpFile() : "";
    
    switch (o.cMode)
    {
        case CalibrateMode::TwoBAM:
        {
            stats.cStats = twoBAM(f1, f2, r1, r2, o);
            calcNorm(r1, r2, stats.cStats, o);
            break;
        }

        case CalibrateMode::Combined:
        {
            stats.cStats = combined(f1, tmp, r2, r3, r4, stats.S1, o); // Write non-sample reads to "tmp"
            calcNorm(r2, r4, stats.cStats, o);
            break;
        }
    }

    stats.tBefore = tBefore(stats.cStats);
    stats.sBefore = sBefore(stats.cStats);

    const auto src = (o.cMode == CalibrateMode::Combined) ? tmp : f2;
    const auto dst = (o.cMode == CalibrateMode::Combined) ? o.work + "/calibrate_combined_calibrated.bam" :
                                                            o.work + "/calibrate_sequin_calibrated.bam";
    
    auto o_ = cloneO(o);
    
    if (o.cMode == CalibrateMode::TwoBAM)
    {
        o_.bam = true;
        o_.flipBefore = true;
    }
    
    // Statistics before normalization
    SCombine(src, stats.S1, o_, &r1); SKallisto(stats.S1, src, "", o_);

    // Build references
    GSplit::buildAF(stats.S1, o); GSplit::buildSL(stats.S1, o);

    // Apply normalization for calibration
    const auto after = o.cMode == CalibrateMode::Combined ?
                          GCalibrate::calibrate(tmp, dst, stats.cStats.norms, r2, r4, stats.cStats.trimmed, stats.cStats, o) :
                          GCalibrate::calibrate(f2,  dst, stats.cStats.norms, r1, r2, stats.cStats.trimmed, stats.cStats, o);
    
    auto o2 = cloneO(o_);
    
    // Statistics after normalization
    SCombine(dst, stats.S2, o2, &r1); SKallisto(stats.S2, dst, "", o2);

    // Build references
    GSplit::buildAF(stats.S2, o); GSplit::buildSL(stats.S2, o);

    // Coverage after calibration
    stats.afterSeqs = o.cMode == CalibrateMode::Combined ?
                          GCalibrate::afterSeqsC(r4, stats.c2v = stats.cStats.c2v, o) :
                          GCalibrate::afterSeqsC(r2, stats.c2v = stats.cStats.c2v, o);
    
    stats.tAfter = tAfter(stats.cStats, after);
    stats.sAfter = sAfter(stats.cStats, after);

    return stats;
}

Coverage GCalibrate::afterSeqsC(const Chr2DInters &r2, std::map<ChrID, std::map<Locus, SampledInfo>> &c2v, const Options &o)
{
    std::vector<Coverage> allAfterSeqsC;
 
    // For each chromosome...
    for (auto &i : c2v)
    {
        const auto &cID = i.first;
        
        // For each region...
        for (auto &j : i.second)
        {
            const auto &l1 = j.first;
         
            std::vector<DInter *> o2;
            
            if (r2.count(cID))
            {
                r2.at(cID).overlap(l1, &o2);
            }
            
            DInter::Stats s1;
            
            if (o2.empty())
            {
                s1.sums = s1.length = s1.mean = s1.p50 = 0;
            }
            else
            {
                std::vector<DInter::Stats> s2;

                for (const auto &j2 : o2)
                {
                    s2.push_back(j2->stats());
                }
                
                assert(!s2.empty());
                
                for (const auto &i : s2)
                {
                    s1.sums += i.sums; s1.aligns += i.aligns; s1.length += i.length; s1.nonZeros += i.nonZeros; s1.zeros += i.zeros;
                    for (const auto &j : i.raws) { s1.raws.push_back(j); }
                }
                
                std::sort(s1.raws.begin(), s1.raws.end());
                
                s1.p50  = SS::quant(s1.raws, 0.50); // Combined statistics
                s1.mean = s1.sums / s1.length;      // Combined statistics
            }

            // Coverage after calibration
            const auto cov = stats2cov(o.meth, s1);
            
            // Alignments after calibration
            const auto aligns = s1.aligns;

            c2v[i.first].at(l1).after  = cov;
            c2v[i.first].at(l1).nAfter = aligns;
            
            // Required for generating summary statistics
            allAfterSeqsC.push_back(cov);
        }
    }

    return SS::mean(allAfterSeqsC);
}

// Total number of alignments before calibration (inside and outside calibration regions)
GCalibrate::GenomeSequins GCalibrate::tBefore(const CalibrateStats &x)
{
    GCalibrate::GenomeSequins r;

    r.nEndo = x.nh;
    r.nSeqs = x.nd;
    
    return r;
}

// Total number of alignments inside calibration regions after calibration
GCalibrate::GenomeSequins GCalibrate::sBefore(const CalibrateStats &x)
{
    GCalibrate::GenomeSequins r;
    
    r.nEndo = x.es.nMap + x.es.nNA;
    r.nSeqs = x.ss.nMap + x.ss.nNA;
    
    return r;
}

// Total number of alignments after calibration (inside and outside calibration regions)
GCalibrate::GenomeSequins GCalibrate::tAfter(const CalibrateStats &x, const ParserBAMBED::Stats &y)
{
    GCalibrate::GenomeSequins r;
    
    // Always equal to before calibration
    r.nEndo = x.nh;
    
    // Anything that has been calibrated and outside calibration region
    r.nSeqs = (y.nNA + y.nMap) + x.nSeqOut;
    
    return r;
}

// Total number of alignments inside calibration regions after calibration
GCalibrate::GenomeSequins GCalibrate::sAfter(const CalibrateStats &x, const ParserBAMBED::Stats &y)
{
    GCalibrate::GenomeSequins r;
    
    r.nEndo = x.es.nMap + x.es.nNA; // Always equal to before calibration
    r.nSeqs = y.nNA + y.nMap;
    
    return r;
}

static void writeCalibrate(const FileName &file, const GCalibrate::Stats &stats, const Options &o)
{
    const auto format = boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8$.2f\t%9$.2f\t%10$.2f\t%11$.2f");
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name"
                                           % "Chrom"
                                           % "Start"
                                           % "End"
                                           % "SampleRead"
                                           % "PreRead"
                                           % "PostRead"
                                           % "SampleCoverage"
                                           % "PreCoverage"
                                           % "PostCoverage"
                                           % "Scale").str());

    // For each chromosome...
    for (const auto &i : stats.c2v)
    {
        // For each variant...
        for (const auto &j : i.second)
        {
            o.writer->write((boost::format(format) % j.second.rID
                                                   % i.first
                                                   % j.first.start
                                                   % j.first.end
                                                   % j.second.nEndo
                                                   % j.second.nBefore
                                                   % j.second.nAfter
                                                   % (o.meth == Method::Percent ? MISSING : toString(j.second.endo))
                                                   % j.second.before
                                                   % j.second.after
                                                   % j.second.norm).str());            
        }
    }
    
    o.writer->close();
}

static void writeReport(const FileName &, const SOptions &o)
{
    if (o.report)
    {
        const auto tmp = o.work + "/" + o.name +  "_features.bed";
        copy(Standard::instance().gen.a1()->src, tmp);
        
        typedef Report::Options::KMOptions KOptions;
        Report::Options o2(o);
        o2.k = std::shared_ptr<KOptions>(new KOptions());
        Report::calibrate(o2);
        removeF(tmp);
    }
}

static void writeSummary(const FileName &file, const FileName &f1, const FileName &f2, const GCalibrate::Stats &stats, const Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto format = "-------SUMMARY STATISTICS\n\n"
                        "-------REFERENCE FILES\n\n"
                        "       Reference variants: %1%\n"
                        "       Reference regions:  %2%\n"
                        "       Reference index:    %3%\n\n"
                        "-------LIBRARY INFORMATION\n\n"
                        "       Version:       %4%\n"
                        "       Instrument ID: %5%\n"
                        "       Run number:    %6%\n"
                        "       Flowcell ID:   %7%\n"
                        "       Lane:          %8%\n\n"
                        "-------USER-SUPPLIED FILES\n\n"
                        "       Input file (first):  %9%\n"
                        "       Input file (second): %10%\n\n"
                        "-------PARAMETERS\n\n"
                        "       K-mer length: %11%\n"
                        "       Threshold:    %12%\n\n"
                        "-------PARTITION SUMMARY\n\n"
                        "       Human reads:             %13% (%14$.4f%%)\n"
                        "       Genome reads:            %15% (%16$.4f%%)\n"
                        "       Ladder reads:            %17% (%18$.4f%%)\n"
                        "       Vector reads:            %19% (%20$.4f%%)\n"
                        "       Structural reads:        %21% (%22$.4f%%)\n"
                        "       Immune reads:            %23% (%24$.4f%%)\n"
                        "       HLA reads:               %25% (%26$.4f%%)\n"
                        "       Information reads:       %27% (%28$.4f%%)\n"
                        "       Mitochondria reads:      %29% (%30$.4f%%)\n"
                        "       Dilution:                %31%%%\n"
                        "       Total:                   %32%\n\n"
                        "-------OUTPUT FASTQ FILES\n\n"
                        "       Human reads path:        %33%\n"
                        "       Genome reads path:       %34%\n"
                        "       Vector reads path:       %35%\n"
                        "       Ladder reads path:       %36%\n"
                        "       Structural reads path:   %37%\n"
                        "       Immune reads path:       %38%\n"
                        "       HLA reads path:          %39%\n"
                        "       Information reads path:  %40%\n"
                        "       Mitochondria reads path: %41%\n\n"
                        "-------Before calibration (within sampling regions)\n\n"
                        "       Sample coverage (average): %42$.2f\n"
                        "       Sequin coverage (average): %43$.2f\n\n"
                        "-------After calibration (within sampling regions)\n\n"
                        "       Sample coverage (average): %44$.2f\n"
                        "       Sequin coverage (average): %45$.2f\n\n"
                        "       Scaling Factor: %46% \u00B1 %47%\n\n"
                        "-------Total alignments (before calibration)\n\n"
                        "       Sample: %48%\n"
                        "       Sequin: %49%\n\n"
                        "-------Total alignments (after calibration)\n\n"
                        "       Sample: %50%\n"
                        "       Sequin: %51%\n\n"
                        "-------Alignments within calibrated regions (before calibration)\n\n"
                        "       Sample: %52%\n"
                        "       Sequin: %53%\n\n"
                        "-------Alignments within calibrated regions (after calibration)\n\n"
                        "       Sample: %54%\n"
                        "       Sequin: %55%\n\n";
    
    //LinearModel af, ld;
    //try { af = stats.af.linear(true, true); } catch(...) {}
    //try { ld = stats.ld.linear(true, true); } catch(...) {}
    
    const auto &r  = Standard::instance().gen;
    const auto &s1 = stats.S1;
    
    #define C(x) (s1.K.c1.count(x) ? s1.K.c1.at(x) : 0)
    #define P(x) (100.0 * s1.K.binP(x))

    const auto gn = C(GR) + C(SO) + C(MS) + C(HP);
    const auto gp = P(GR) + P(SO) + P(MS) + P(HP);
    assert(gp >= 0 && gp <= 100);
    
    extern FASTQ __KFQInfo__;
    const auto f = __KFQInfo__.format();
    
    o.writer->write((boost::format(format) % r.v1()->src // 1
                                           % (o.bam ? r.r1()->src + " and " + r.r2()->src : r.r1()->src)
                                           % o.index
                                           % SVersion(Standard::instance().gen, s1.K) // 4
                                           % __KFQInfo__.inst(f)                      // 5
                                           % __KFQInfo__.run(f)                       // 6
                                           % __KFQInfo__.flow(f)                      // 7
                                           % __KFQInfo__.lane(f)                      // 8
                                           % f1                 // 9
                                           % f2                 // 10
                                           % o.k                // 11
                                           % o.rule             // 12
                                           % S0(C(ES))          // 13
                                           % S2(P(ES))          // 14
                                           % S0(gn)             // 15
                                           % S2(gp)             // 16
                                           % S0(C(LD))          // 17
                                           % S2(P(LD))          // 18
                                           % S0(C(VC))          // 19
                                           % S2(P(VC))          // 20
                                           % S0(C(SV))          // 21
                                           % S2(P(SV))          // 22
                                           % S0(C(IM))          // 23
                                           % S2(P(IM))          // 24
                                           % S0(C(HL))          // 25
                                           % S2(P(HL))          // 26
                                           % S0(C(IF))          // 27
                                           % S2(P(IF))          // 28
                                           % S0(C(MT))          // 29
                                           % S2(P(MT))          // 30
                                           % (100.0 * s1.dil()) // 31
                                           % s1.K.total()       // 32
                                           % (o.work + "/split_sample*") // 33
                                           % (o.work + "/split_sequin*") // 34
                                           % (o.work + "/split_vector*") // 35
                                           % (o.work + "/split_ladder*") // 36
                                           % (o.work + "/split_sv*")     // 37
                                           % (o.work + "/split_immune*") // 38
                                           % (o.work + "/split_hla*")    // 39
                                           % (o.work + "/split_info*")   // 40
                                           % (o.work + "/split_mito*")   // 41
                                           % (o.meth == Method::Percent ? "-" : toString(stats.cStats.meanBEndo())) // 42
                                           % stats.cStats.meanBSeqs() // 43
                                           % (o.meth == Method::Percent ? "-" : toString(stats.cStats.meanBEndo()))
                                           % stats.afterSeqs          // 45
                                           % stats.cStats.normMean()  // 46
                                           % stats.cStats.normSD()    // 47
                                           % stats.tBefore.nEndo      // 48
                                           % stats.tBefore.nSeqs      // 49
                                           % (o.meth == Method::Percent ? "-" : toString(stats.tAfter.nEndo))
                                           % stats.tAfter.nSeqs       // 51
                                           % (o.meth == Method::Percent ? "-" : toString(stats.sBefore.nEndo))
                                           % stats.sBefore.nSeqs      // 53
                                           % (o.meth == Method::Percent ? "-" : toString(stats.sAfter.nEndo))
                                           % stats.sAfter.nSeqs       // 55
                                    ).str());
    o.writer->close();
}

static void commonReport(const GCalibrate::Stats &stats, const GCalibrate::Options &o)
{
    // Generating calibrate_calibrate.tsv
    writeCalibrate("calibrate_calibrate.tsv", stats, o);

    // Generating calibrate_sequin.tsv
    GSplit::writeQuin("calibrate_sequin.tsv", stats.S1, o);

    // Generating calibrate_sequin.tsv
    GSplit::writeQuin("calibrate_calibrated_sequin.tsv", stats.S2, o);
    
    // Generating calibrate_somatic.R
    writeSomatic("calibrate_somatic.R", "calibrate_calibrated_sequin.tsv", o);
    
    // Generating calibrate_reads.tsv
    SWriteReads(Product::Genomics, "calibrate_read.tsv", stats.S2, o);
    
    // Synthetic calibration
    SCalibSynthetic(stats.S1, GSplit::analyze, o.work, o, Standard::instance().gen.l3());

    // Generating calibrate_report.html
    writeReport("calibrate_report.html", o);    
}

void GCalibrate::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, "", o);
    
    // Generating calibrate_summary.stats
    writeSummary("calibrate_summary.stats", file, "", stats, o);
    
    // Everything else will be the same
    commonReport(stats, o);
}

void GCalibrate::report(const FileName &hr, const FileName &dr, const Options &o)
{
    const auto stats = analyze(hr, dr, o);
    
    // Generating calibrate_summary.stats
    writeSummary("calibrate_summary.stats", hr, dr, stats, o);

    // Everything else will be the same
    commonReport(stats, o);
}
