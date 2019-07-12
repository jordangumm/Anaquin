#include "writers/vcf_writer.hpp"
#include "Genomics/g_variant.hpp"

using namespace Anaquin;

enum Mode
{
    HumanVCF = 1,
    DecoyVCF = 2,
    HumanAndDecoy = 3
};

const auto muts = std::set<Variation>
{
    Variation::SNP,
    Variation::Deletion,
    Variation::Insertion,
};

inline bool isShort(const SequinID &x)
{
    return isSubstr(x, "GV_") ||
           isSubstr(x, "DV_") ||
           isSubstr(x, "HP_") ||
           isSubstr(x, "CV_") ||
           isSubstr(x, "CM_") ||
           isSubstr(x, "CL_") ||
           isSubstr(x, "MS_");
}

static void analyzeMode(Mode mode, const GVariant &v, GVariant::Stats &stats, const FileName &file, const GVariant::Options &o, std::shared_ptr<VCFLadder> vl)
{
    auto w1 = (HumanVCF & mode) ? std::shared_ptr<VCFWriter>(new VCFWriter()) : nullptr;
    auto wT = (DecoyVCF & mode) ? std::shared_ptr<VCFWriter>(new VCFWriter()) : nullptr;
    auto wF = (DecoyVCF & mode) ? std::shared_ptr<VCFWriter>(new VCFWriter()) : nullptr;

    if (w1) { w1->open(o.work + "/genome_files/" + o.name + "_sample.vcf"); }
    if (wT) { wT->open(o.work + "/genome_files/" + o.name + "_TP.vcf");      }
    if (wF) { wF->open(o.work + "/genome_files/" + o.name + "_FP.vcf");      }

    const auto &r = Standard::instance().gen;

    const auto r2 = r.r2()->inters(); // Human regions including edges
    const auto r3 = r.r3()->inters(); // Decoy regions without edges
    const auto r4 = r.r4()->inters(); // Decoy regions including edges

    // We'll need it for false negatives
    std::set<long> keys;

    o.analyze(file);
    
    ParserVCF::parse(file, [&](const Variant &x)
    {
        DInter *d = nullptr;

        if (HumanVCF & mode)
        {
            assert(w1);
            
            if ((d = contains(r2, x.cID, x.l)) && isShort(d->name()))
            {
                w1->write(x.hdr, x.line);
                auto t = x;
                t.name = d->name();
                stats.hs.vs.insert(t);
            }
        }
        
        if (!d && (DecoyVCF & mode) && (d = contains(r4, x.cID, x.l)))
        {
            auto findMatch = [&](const Variant &q)
            {
                VCFMatch m;
                
                m.qry = q;
                m.var = nullptr;
                
                // Can we match by position?
                if ((m.var = r.v1()->data.findVar(q.cID, q.l)))
                {
                    // Match by reference allele?
                    m.ref = m.var->ref == q.ref;
                    
                    // Match by alternative allele?
                    m.alt = m.var->alt == q.alt;
                }
                
                if (m.var)
                {
                    m.rID = m.var->name;
                    A_ASSERT(!m.rID.empty());
                }
                else
                {
                    GIntervals<> inters;
                    
                    try
                    {
                        // Search where the FPs are (finding shoudl ignore edges)
                        inters = r.r3()->ginters(x.cID);
                        
                        const auto m2 = inters.contains(x.l);
                        
                        // Can we find the corresponding region for the FP?
                        if (m2)
                        {
                            m.rID = m2->name();
                            A_ASSERT(!m.rID.empty());
                        }
                    }
                    catch (...) {}
                }
                
                return m;
            };
            
            auto m = findMatch(x);
            
            if (!m.rID.empty() && !v.isValid(m.rID))
            {
                return;
            }
            
            // FP because not in germline regions?
            const auto fFP = m.var && !v.isValid(m.var->name);
            
            if (fFP)
            {
                // Let's make it unmatched
                m.var = nullptr;
            }
            
            // Matched if both the position and alleles agree
            const auto matched = m.var && m.ref && m.alt;
            
            if (matched)
            {
                keys.insert(x.key());
                wT->write(x.hdr, x.line);
                stats.ds.tps.push_back(m);
                A_ASSERT(!std::isnan(r.af(m.var->name)));
            }
            else
            {
                if (!fFP)
                {
                    wF->write(x.hdr, x.line);
                    stats.ds.fps.push_back(m);
                }
            }
        }
    });
    
    if (DecoyVCF & mode)
    {
        /*
         * Determine the classification performance
         */
        
        auto forTP = [&]()
        {
            for (auto i = 0u; i < stats.ds.tps.size(); i++)
            {
                // Overall performance
                stats.ds.oc.tp()++;
            }
        };
        
        auto forFP = [&]()
        {
            for (auto i = 0u; i < stats.ds.fps.size(); i++)
            {
                // Overall performance
                stats.ds.oc.fp()++;
            }
        };
        
        forTP();
        forFP();
        
        for (auto &mut : muts)
        {
            stats.ds.oc.nr() += r.nType(vl, mut);
        }
        
        stats.ds.oc.fn() = stats.ds.oc.nr() - stats.ds.oc.tp();
        A_ASSERT(stats.ds.oc.nr() >= stats.ds.oc.fn());
        
        /*
         * Work out the false negatives
         */
        
        for (const auto &i : r.v1()->data.vars())
        {
            if (v.isValid(i.name) && !stats.ds.findTP(i.name))
            {
                VCFMatch m;
                
                m.var = r.v1()->data.findVar(i.cID, i.l);
                m.rID = i.name;
                A_ASSERT(m.var);
                
                if (!keys.count(m.var->key()))
                {
                    stats.ds.fns.push_back(m);
                }
            }
        }
    }
}

GVariant::Stats GVariant::analyze(Stats &stats, const GVariant &x, const FileName &f1, const FileName &f2, const Options &o, std::shared_ptr<VCFLadder> vl)
{
    o.info("Edge: " + toString(o.edge));
    o.info("Build: " + toString(o.build));
    o.info("Variant: " + Standard::instance().gen.v1()->src);

    if (o.combined)
    {
        assert(f2.empty());
        analyzeMode(HumanAndDecoy, x, stats, f1, o, vl);
    }
    else
    {
        if (!f1.empty()) { analyzeMode(HumanVCF, x, stats, f1, o, vl); }
        analyzeMode(DecoyVCF, x, stats, f2, o, vl);
    }

    o.info("TP: " + toString(stats.ds.oc.tp()));
    o.info("FP: " + toString(stats.ds.oc.fp()));
    o.info("FN: " + toString(stats.ds.fns.size()));

    return stats;
}
