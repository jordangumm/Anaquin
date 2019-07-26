#include <fstream>
#include "tools/random.hpp"
#include "stats/ss/stats.hpp"
#include "tools/calibrate.hpp"
#include "Genomics/Genomics.hpp"
#include "parsers/parser_fq.hpp"
#include "writers/bam_writer.hpp"
#include "Genomics/g_synthetic.hpp"

using namespace Anaquin;

static std::map<StandardID, KMCoverage> coverage(const KStats &ks)
{
    std::map<StandardID, std::vector<KMCoverage>> x1;
    std::map<StandardID, KMCoverage> x2;
    
    // Calibrate by number of reads
    for (auto &i : ks.sqc)
    {
        if (isSubstr(i.first, "LD_"))
        {
            x1[GSeq2Std(i.first)].push_back(i.second);
        }
    }

    for (auto &i : x1)
    {
        const auto sum = SS::Internal::sum(i.second);

        if (sum)
        {
            x2[i.first] = sum;
        }
    }

    return x2;
}

static Selection selectPool(Proportion p, const std::map<StandardID, KMCoverage> &x, const WriterOptions &o, LadderInternal::Results &r)
{
    assert(!std::isnan(p) && p >= 0.0 && p <= 1.0);
    Selection rnd;

    // Pool size for all ladders
    const auto pool = sum(x);
    
    // Proportion of the pool
    r.targetC = p * pool;
    
    std::vector<double> scales;
    
    for (const auto &i : x)
    {
        // Discard if less than the target
        const auto scale = i.second < r.targetC ? 1.0 : 1.0 - (r.targetC / i.second);
        
        scales.push_back(scale);
        
        // Scale each individual ladder to the target
        rnd[i.first] = std::shared_ptr<RandomSelection>(new RandomSelection(scale));
    }

    // Aveage calibration scaling factor
    r.meanS = SS::mean(scales);

    return rnd;
}

LadderInternal::Results LadderInternal::calibrate(const KStats &stats,
                                          const Path &src,
                                          const KOptions &ko,
                                          Proportion p,
                                          const WriterOptions &o,
                                          const Label &ll)
{
    LadderInternal::Results r;
    r.cov = coverage(stats);

    const auto i1 = !ko.writeBAM() ? (src + "/" + ll + "_ladder_1.fq.gz") : (src + "/" + ll + "_ladder.bam");
    const auto i2 = !ko.writeBAM() ? (src + "/" + ll + "_ladder_2.fq.gz") : "";

    if (r.cov.empty())
    {
        o.warn("Ladder reads not found. Nothing to calibrate.");
        r.r1 = i1; r.r2 = i2; // No calibration
        return r;
    }
    else if (p == NO_CALIBRATION)
    {
        o.info("Second stage ladder calibration skipped");
        r.r1 = i1; r.r2 = i2; // No calibration
        return r;
    }
    
    // Calibration by pool can only be a proportion
    assert(p >= 0.0 && p <= 1.0);
    
    for (const auto &i : r.cov)
    {
        o.logInfo(i.first + ": " + toString(i.second));
    }

    auto calib = (ko.writeBAM()) ? Calibrator::createBAM(i1, o.work + "/" + ll + "_ladder_calibrated.bam") :
                                   Calibrator::createFQ(i1, i2, o.work + "/" + ll + "_ladder_calibrated_1.fq.gz", o.work + "/" + ll + "_ladder_calibrated_2.fq.gz");

    // Normalizate against the all ladder reads (pooled)
    auto rnd = selectPool(p, r.cov, o, r);
    
    auto c = calib->calibrate(stats, rnd, o);
    r.r1 = c.o1; r.r2 = c.o2; r.after = c.n;
    return r;
}
