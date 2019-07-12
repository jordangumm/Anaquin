#include <fstream>
#include "ss/stats.hpp"
#include "tools/random.hpp"
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

static Selection scaling(const std::map<StandardID, KMCoverage> &x_, const WriterOptions &o, GSynthetic::Results &r)
{
    Selection rnd;
    
    auto x = x_;
    const auto it = std::min_element(x.begin(), x.end(),
                                    [](const decltype(x)::value_type &l, const decltype(x)::value_type &r) -> bool { return l.second < r.second; });

    r.target  = it->first;
    r.targetC = it->second;
    
    for (auto &i : x)
    {
        assert(i.second >= r.targetC);
        rnd[i.first] = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - (r.targetC / i.second)));
    }

    return rnd;
}

GSynthetic::Results GSynthetic::calibrate(const KStats &stats,
                                          const Path &src,
                                          const KOptions &ko,
                                          const WriterOptions &wo,
                                          const Label &lab)
{
    GSynthetic::Results r;
    r.cov = coverage(stats);
    
    const auto i1 = !ko.writeBAM() ? (src + "/" + lab + "_ladder_1.fq.gz") : (src + "/" + lab + "_ladder.bam");
    const auto i2 = !ko.writeBAM() ? (src + "/" + lab + "_ladder_2.fq.gz") : "";

    if (r.cov.empty())
    {
        wo.warn("Ladder reads not found. Nothing to calibrate");
        r.r1 = i1; r.r2 = i2; // No calibration
        return r;
    }
    
    wo.logInfo("Running synthetic calibration");
    
    // Random generator
    auto rnd = scaling(r.cov, wo, r);

    for (const auto &i : r.cov)
    {
        wo.logInfo(i.first + ": " + toString(i.second));
    }
    
    auto calib = (ko.writeBAM()) ? Calibrator::createBAM(i1,
                                                         wo.work + "/" + lab + "_ladder_calibrated.bam") :
                                   Calibrator::createFQ (i1, i2,
                                                         wo.work + "/" + lab + "_ladder_calibrated_1.fq.gz",
                                                         wo.work + "/" + lab + "_ladder_calibrated_2.fq.gz");

    auto c = calib->calibrate(stats, rnd, wo);
    r.r1 = c.o1; r.r2 = c.o2;
    return r;
}
