#ifndef G_VARIANT_HPP
#define G_VARIANT_HPP

#include "stats/analyzer.hpp"
#include "Genomics/Genomics.hpp"
#include "writers/vcf_writer.hpp"

namespace Anaquin
{
    struct VCFMatch
    {
        // The called variant
        Variant qry;
        
        // Sequin matched by position?
        const Variant *var = nullptr;
        
        // Matched by variant allele? Only if position is matched.
        bool alt;
        
        // Matched by reference allele? Only if position is matched.
        bool ref;
        
        // Does the variant fall into one of the reference regions?
        SequinID rID;
    };

    template <typename T, typename O> void writeFN(const FileName &file, const T &g, const O &o)
    {
        o.generate(file);
        
        auto wFN = VCFWriter();
        wFN.open(o.work + "/genome_files/" + file);
        
        std::set<long> keys;
        for (const auto &i : g) { keys.insert(i.var->key()); }
        
        ParserVCF::parse(Standard::instance().gen.v1()->src, [&](const Variant &x)
        {
            if (keys.count(x.key()))
            {
                wFN.write(x.hdr, x.line);
            }
        });
    }

    struct GVariant
    {
        struct HStats
        {
            std::set<Variant> vs;
            
            /*
             * Caller specific fields
             */
            
            std::map<std::string, std::map<long, int>>   si;
            std::map<std::string, std::map<long, float>> sf;
        };
        
        struct DStats
        {
            std::vector<VCFMatch> tps, fns, fps;
            
            // Overall performance
            Confusion oc;
            
            inline const VCFMatch * findTP(const SequinID &id) const
            {
                for (auto &i : tps)
                {
                    if (i.var->name == id)
                    {
                        return &i;
                    }
                }
                
                return nullptr;
            }
            
            /*
             * Caller specific fields
             */
            
            std::map<std::string, std::map<long, int>>   si;
            std::map<std::string, std::map<long, float>> sf;
        };
        
        struct Stats
        {
            HStats hs;
            DStats ds;
        };
        
        struct Options : public AnalyzerOptions<HumanAssembly>
        {
            Base edge;
            
            FileName uBED;
          
            bool combined;
            
            // Generating HTML report?
            bool report;
        };
        
        virtual bool isValid(const SequinID &) const = 0;
        
        static Stats analyze(Stats &, const GVariant &, const FileName &, const FileName &, const Options &, std::shared_ptr<VCFLadder>);
    };
}

#endif
