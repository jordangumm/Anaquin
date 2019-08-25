#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include "tools/report.hpp"
#include "tools/random.hpp"
#include "parsers/parser_fa.hpp"
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

inline void addInsert(ParserBAM::Data &x, SAlignStats &align, bool isSeq)
{
    if (x.isPrimary)
    {
        x.lMateID(); x.lMatePos();
        
        if (x.cID == x.mID)
        {
            // Starting position of this read
            auto x1 = x.l.start;
            
            // Starting position of mate
            auto x2 = x.mPos;
            
            if (x2 < x1)
            {
                // Always make sure sure x1 before x2
                std::swap(x1, x2);
            }
            
            // Insertion size
            const auto ins = x2 - x1;

            // Update frequency table
            align.ins[isSeq ? GR : ES][ins]++;
        }
    }
}

static ParserBAMBED::Stats readSamp(const FileName &file, SAlignStats &align, const Chr2DInters &r, GCalibrate::CalibrateStats &stats, const GCalibrate::Options &o)
{
    return ParserBAMBED::parse(file, r, [&](ParserBAM::Data &x, const ParserBAM::Info &info, const DInter *m)
    {
        if (info.p && !(info.p % 1000000)) { o.wait(toString(info.p)); }
        
        if (x.mapped)
        {
            stats.nh++;
        }
        
        // Update insertion size for sample but only within regions
        if (x.mapped && m)
        {
            addInsert(x, align, false);
        }
        
        return ParserBAMBED::Response::OK;
    });
}

static ParserBAMBED::Stats readSeq(const FileName &file, const Chr2DInters &r, GCalibrate::CalibrateStats &stats, const GCalibrate::Options &o)
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
                
                g1.p50  = quant(g1.raws, 0.50); // Combined statistics
                s1.p50  = quant(s1.raws, 0.50); // Combined statistics
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
                                                SAlignStats &align,
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
        if (i.p && !(i.p % 100000)) { o.wait(toString(i.p)); }

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
        
        if (!decoy)
        {
            // Update insertion size for sample
            addInsert(x, align, false);
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
                                              SAlignStats &align,
                                              const Options &o)
{
    GCalibrate::CalibrateStats stats;
    
    // Check sample alignments before calibration
    stats.es = readSamp(endo, align, r2, stats, o);
    
    // Calculate trimming on sequin alignments
    trimming(seqs, r1, stats, o);
    
    // Check sequin alignments before calibration
    stats.ss = readSeq(seqs, r2, stats, o);
    
    return stats;
}

struct PicardData
{
    // Chromsome in FASTA
    std::map<ChrID, Sequence> chrs;
    
    /*
     * Mapping from BED region names to FASTA names (chrs). Important because they might not be fully compatible.
     */
    
    std::map<ChrID, ChrID> m;
};

static PicardData initPicard(bool isDecoy, const std::map<ChrID, DIntervals<>> &r1, const GCalibrate::GCalibrate::Options &o)
{
    /*
     * Sequin names in FASTA and BED can be inconsistent. Furthermore, we would like to construct mapping for only reference
     * sequins ("_R"). Our Picard won't work on "_A".
     */
    
    PicardData r;
    
    if (o.cMode == GCalibrate::CalibrateMode::TwoBAM)
    {
        ParserFA::parse(Reader(o.decoy), [&](const ParserFA::Data &x) {
            if (x.id == GENOMICS_DECOY_CHROM)
            {
                r.chrs[x.id] = x.seq;
            }
        });
        
        // Only chrQS is supported
        assert(r.chrs.size() == 1);
    }

    return r;
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

    // Initalize data requied for Picard
    const auto pd = initPicard(true, r1, o);
    
    // Picard after calibration
    //stats.P2 = o.cMode == CalibrateMode::Combined ? std::shared_ptr<Picard>(new Picard(pd.chrs)) : nullptr;
    
    // BAM file for combined mode
    const auto tmp = o.cMode == CalibrateMode::Combined ? tmpFile() : "";
    
    // Only used for sample
    auto align = SAlignStats(); align.ins[ES];
    
    switch (o.cMode)
    {
        case CalibrateMode::TwoBAM:
        {
            stats.cStats = twoBAM(f1, f2, r1, r2, align, o);
            calcNorm(r1, r2, stats.cStats, o);
            break;
        }

        case CalibrateMode::Combined:
        {
            stats.cStats = combined(f1, tmp, r2, r3, r4, stats.S1, align, o); // Write non-sample reads to "tmp"
            calcNorm(r2, r4, stats.cStats, o);
            break;
        }
    }
    
    // Only used for sample reads (sequin reads below)
    assert(align.ins[GR].empty());

    stats.tBefore = tBefore(stats.cStats);
    stats.sBefore = sBefore(stats.cStats);

    const auto src = (o.cMode == CalibrateMode::Combined) ? tmp : f2;
    const auto dst = o.work + "/calibrate_sequins.bam";
    
    auto o_ = cloneO(o);
    
    if (o.cMode == CalibrateMode::TwoBAM)
    {
        o_.bam = true;
        o_.flipBefore = true;
    }
    
    const auto rc = false; // TODO: o.cMode == CalibrateMode::Combined;
    
    // Statistics before normalization
    stats.A1 = SAlign(src, stats.S1, o_, &r1, rc); SKallisto(stats.S1, src, "", o_);

    // There shouldn't be any sample read
    stats.A1.ins[ES]; assert(stats.A1.ins[ES].empty()); stats.A1.ins[ES] = align.ins[ES];

    // Build references
    GSplit::buildAF(stats.S1, o); GSplit::buildSL(stats.S1, o);

    // Apply normalization for calibration
    const auto after = o.cMode == CalibrateMode::Combined ?
                          GCalibrate::calibrate(tmp, dst, stats.cStats.norms, r2, r4, stats.cStats.trimmed, stats.cStats, o) :
                          GCalibrate::calibrate(f2,  dst, stats.cStats.norms, r1, r2, stats.cStats.trimmed, stats.cStats, o);
    
    auto o2 = cloneO(o_);
    
    // Statistics after normalization
    stats.A2 = SAlign(dst, stats.S2, o2, &r1, rc); SKallisto(stats.S2, dst, "", o2);

    // There shouldn't be any sample reads
    stats.A2.ins[ES]; assert(stats.A2.ins[ES].empty()); stats.A2.ins[ES] = align.ins[ES];

    /*
     * Work out sequecing profile from the calibrated BAM (not before calibration)
     */

    PicardOption op;
    op.ignoreSkipClip = true; // Skips and clipping not supported in this release
    
/*
    if (o.cMode == CalibrateMode::Combined)
    {
        ParserBAM::parse(dst, [&](ParserBAM::Data &x, const ParserBAM::Info &)
        {
            if (x.cID == GENOMICS_DECOY_CHROM)
            {
                x.lCigar(); x.lSeq(); stats.P2->analyze(GENOMICS_DECOY_CHROM, x, op);
            }
        });
    }
*/
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
                
                s1.p50  = quant(s1.raws, 0.50); // Combined statistics
                s1.mean = s1.sums / s1.length;  // Combined statistics
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
    o.writer->write((boost::format(format) % "NAME"
                                           % "CHROM"
                                           % "START"
                                           % "END"
                                           % "SAMPLE_READ"
                                           % "PRE_READ"
                                           % "POST_READ"
                                           % "SAMPLE_COVERAGE"
                                           % "PRE_COVERAGE"
                                           % "POST_COVERAGE"
                                           % "SCALE").str());

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
        // R-script for insertion size
        writeInsertR("report_files/calibrate_insert.R", "calibrate_insert.tsv", o);

        // Generating calibrate_somatic.R
        writeSomatic("report_files/calibrate_somatic.R", "calibrate_sequin.tsv", o);

        const auto tmp = o.work + "/" + o.name +  "_features.bed";
        copy(Standard::instance().gen.a1()->src, tmp);
        
        typedef Report::Options::KMOptions KOptions;
        Report::Options o2(o);
        o2.k = std::shared_ptr<KOptions>(new KOptions());
        Report::calibrate(o2);
        removeF(tmp);
    }
}

static void writePicard(const Picard &p, const FileName &file, const Options &o)
{
    if (o.cMode == GCalibrate::CalibrateMode::Combined)
    {
        o.generate(file);
        o.writer->open(file);
        o.writer->write(p.report());
        o.writer->close();
    }
}

static void writeSummary(const FileName &file, const FileName &f1, const FileName &f2, const GCalibrate::Stats &stats, const Options &o)
{
    const auto tsv = o.isSCalib() ? o.work + "/calibrate_sequin_calibrated.tsv" : o.work + "/calibrate_sequin.tsv";
    
    // Generating sequin abundance table
    writeSTable(tsv, "calibrate_sequin_table.tsv", o, 6, 6, 6, "EXP_FREQ", "OBS_FREQ");

    const auto tmp = tmpFile();
    RGrep(tsv, tmp, "LABEL", "Somatic"); const auto l2 = RLinear(tmp, "NAME", "EXP_FREQ", "OBS_FREQ").linear();

    const auto f = "SEQUIN REPORT:                 %1%\n\n"
                   "REFERENCE FILES\n"
                   "Reference index:               %2%\n\n"
                   "LIBRARY INFORMATION\n"
                   "Version:                       %3%\n"
                   "Instrument ID:                 %4%\n"
                   "Run number:                    %5%\n"
                   "Flowcell ID:                   %6%\n"
                   "Lane:                          %7%\n\n"
                   "USER-SUPPLIED FILES\n"
                   "Input file (first):            %8%\n"
                   "Input file (second):           %9%\n\n"
                   "ANAQUIN PARAMETERS\n"
                   "K-mer length:                  %10%\n"
                   "Threshold:                     %11%\n\n"
                   "PARTITION SUMMARY (BEFORE CALIBRATION)\n"
                   "Sample reads:                  %12% (%13%%%)\n"
                   "Sequin reads:                  %14% (%15%%%)\n"
                   "Ladder reads:                  %16% (%17%%%)\n"
                   "Structural reads:              %18% (%19%%%)\n"
                   "Immune reads:                  %20% (%21%%%)\n"
                   "HLA reads:                     %22% (%23%%%)\n"
                   "Mitochondria reads:            %24% (%25%%%)\n"
                   "Vector reads:                  %26% (%27%%%)\n"
                   "Information reads:             %28% (%29%%%)\n"
                   "Dilution:                      %30%%%\n"
                   "Total reads:                   %31%\n\n"
                   "OUTPUT FILES\n"
                   "Sample reads path:             %32%\n"
                   "Sequin reads path:             %33%\n"
                   "Ladder reads path:             %34%\n"
                   "Structural reads path:         %35%\n"
                   "Immune reads path:             %36%\n"
                   "Mitochondria reads path:       %37%\n"
                   "HLA reads path:                %38%\n"
                   "Information reads path:        %39%\n"
                   "Vector reads path:             %40%\n\n"
                   "CALIBRATION SUMMARY\n"
                   "Before calibration (within sampling regions)\n"
                   "Sample alignments:             %41%\n"
                   "Sequin alignments:             %42%\n"
                   "Sample coverage (mean):        %43%\n"
                   "Sequin coverage (mean):        %44%\n\n"
                   "After calibration (within sampling regions)\n"
                   "Scaling Factor:                %45% \u00B1 %46%\n"
                   "Sample alignments:             %47%\n"
                   "Sequin alignments:             %48%\n"
                   "Sample coverage (mean):        %49%\n"
                   "Sequin coverage (mean):        %50%\n\n"
    /*
                   "SEQUENCING PROFILE:\n"
                   "Mismatch total:                %51%\n"
                   "Mismatch per KB:               %52%\n"
                   "Insertions total:              %53%\n"
                   "Insertions per KB:             %54%\n"
                   "Deletions total:               %55%\n"
                   "Deletions per KB:              %56%\n"
                   "Sequencing error table:        %57%\n\n"
                   "Mean read insert size (nt):    %58%\n"
                   "Insert size table:             %59%\n\n"
     */
                   "SEQUIN SOMATIC - LIBRARY QUALITY\n"
                   "Slope:                         %51%\n"
                   "R2:                            %52%\n"
                   "Sequin table:                  %53%";

    const auto &r  = Standard::instance().gen;
    const auto &s1 = stats.S1;
    
    #define C(x) (s1.K.c1.count(x) ? s1.K.c1.at(x) : 0)
    #define P(x) (100.0 * s1.K.binP(x))

    const auto gn = C(GR) + C(SO) + C(MS) + C(HP);
    const auto gp = P(GR) + P(SO) + P(MS) + P(HP);
    assert(gp >= 0 && gp <= 100);
    
    extern FASTQ __KFQInfo__;
    const auto fo = __KFQInfo__.format();
/*
    // Only chrQS is supported
    assert(!stats.P2 || stats.P2->seqs.size() == 1);
    
    // Total number of mismatches (valid only if combined mode)
    const auto N = stats.P2 ? (stats.P2->totalS() - stats.P2->sumSB(SNPBin::Match)) : 0;
    
    // Total length (valid only if combined mode)
    const auto L = stats.P2 ? stats.P2->seqs.at(GENOMICS_DECOY_CHROM).size() : 0;

    // Total number of mismatches per KB (valid only if combined mode)
    const auto NKB = stats.P2 ? ((N / (double) L) * 1000.0) : 0.0;
    
    const auto D = stats.P2 ? stats.P2->sumDL() : 0;
    const auto I = stats.P2 ? stats.P2->sumIL() : 0;
    
    const auto DL = stats.P2 ? ((D / (double) L) * 1000.0) : 0.0;
    const auto IL = stats.P2 ? ((I / (double) L) * 1000.0) : 0.0;
*/
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % date()                  // 1
                                      % o.index                 // 2
                                      % SVersion(r, stats.S1.K) // 3
                                      % __KFQInfo__.inst(fo)    // 4
                                      % __KFQInfo__.run(fo)     // 5
                                      % __KFQInfo__.flow(fo)    // 6
                                      % __KFQInfo__.lane(fo)    // 7
                                      % f1                      // 8
                                      % (f2.empty() ? "-" : f2) // 9
                                      % o.k                     // 10
                                      % o.rule                  // 11
                                      % C(ES)                   // 12
                                      % S2(100.0 * stats.S1.K.binP(ES))
                                      % C(GR)                // 14
                                      % S2(100.0 * stats.S1.K.binP(GR))
                                      % C(LD)                // 16
                                      % S2(100.0 * stats.S1.K.binP(LD))
                                      % C(SV)               // 18
                                      % S2(100.0 * stats.S1.K.binP(SV))
                                      % C(IM)               // 20
                                      % S2(100.0 * stats.S1.K.binP(IM))
                                      % C(HL)               // 22
                                      % S2(100.0 * stats.S1.K.binP(HL))
                                      % C(MI)               // 24
                                      % S2(100.0 * stats.S1.K.binP(MI))
                                      % C(VC)               // 26
                                      % S2(100.0 * stats.S1.K.binP(VC))
                                      % C(IF)               // 28
                                      % S2(100.0 * stats.S1.K.binP(IF))
                                      % S2(100.0 * s1.dil()) // 30
                                      % stats.S1.K.total()   // 31
                                      % (o.work + "/calibrate_sample*") // 32
                                      % (o.work + "/calibrate_sequin*") // 33
                                      % (o.work + "/calibrate_ladder*") // 34
                                      % (o.work + "/calibrate_sv*")     // 35
                                      % (o.work + "/calibrate_immune*") // 36
                                      % (o.work + "/calibrate_mito*")   // 37
                                      % (o.work + "/calibrate_hla*")    // 38
                                      % (o.work + "/calibrate_info*")   // 39
                                      % (o.work + "/calibrate_vector*") // 40
                                      % stats.tBefore.nEndo         // 41
                                      % stats.tBefore.nSeqs         // 42
                                      % (o.meth == Method::Percent ? MISSING : toString(stats.cStats.meanBEndo())) // 43
                                      % stats.cStats.meanBSeqs()    // 44
                                      % S2(stats.cStats.normMean()) // 45
                                      % S2(stats.cStats.normSD())   // 46
                                      % (o.meth == Method::Percent ? MISSING : toString(stats.tAfter.nEndo))
                                      % stats.tAfter.nSeqs          // 48
                                      % (o.meth == Method::Percent ? MISSING : S2(stats.cStats.meanBEndo()))
                                      % S2(stats.afterSeqs)         // 50
                                      //% (stats.P2 ? std::to_string(N) : MISSING)
                                      //% (stats.P2 ? std::to_string(NKB) : MISSING)
                                      //% D  // 53
                                      //% DL // 54
                                      //% I  // 55
                                      //% IL // 56
                                      //% (o.work +  "/calibrate_seq_errors.tsv") // 57
                                      //% (stats.A2.ins.count(Bin::GR) ? toString(stats.A2.mInsert()) : MISSING)
                                      //% (o.work +  "/calibrate_insert.tsv") // 59
                                      % replaceNA(l2.m)  // 50
                                      % replaceNA(l2.R2) // 51
                                      % (o.work + "/calibrate_sequin_table.tsv") // 52
                     ).str());
    o.writer->close();
}

static void commonReport(const GCalibrate::Stats &stats, const GCalibrate::Options &o)
{
    // Generating calibrate_regions.tsv
    writeCalibrate("calibrate_regions.tsv", stats, o);

    // Generating calibrate_reads.tsv
    SWriteReads(Product::Genomics, "calibrate_read.tsv", stats.S2, o);
    
    // Synthetic calibration
    StageTwoLadderCalibration(stats.S1, GSplit::analyze, o.work, o, Standard::instance().gen.l3());

    // TSV for insertion size
    writeInsertTSV("calibrate_insert.tsv", stats.A2, o);
    
    // Generate sequencing profile
    //writePicard(*stats.P2.get(), "calibrate_seq_errors.tsv", o);

    // Generating HTML report
    writeReport("calibrate_report.html", o);    
}

void GCalibrate::report(const FileName &file, const Options &o)
{
    const auto stats = analyze(file, "", o);
    
    // Generating performance after calibration (stats.S1 for before calibration and not written)
    GSplit::writeQuin("calibrate_sequin.tsv", stats.S2, o);

    // Generating calibrate_summary.txt
    writeSummary("calibrate_summary.txt", file, "", stats, o);
    
    // Everything else will be the same
    commonReport(stats, o);
}

void GCalibrate::report(const FileName &hr, const FileName &dr, const Options &o)
{
    const auto stats = analyze(hr, dr, o);
    
    // Generating performance after calibration (stats.S1 for before calibration)
    GSplit::writeQuin("calibrate_sequin.tsv", stats.S2, o);

    // Generating calibrate_summary.txt
    writeSummary("calibrate_summary.txt", hr, dr, stats, o);

    // Everything else will be the same
    commonReport(stats, o);
}
