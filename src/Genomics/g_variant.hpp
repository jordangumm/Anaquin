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
            
            // Required for common operations
            Label base;
        };
        
        virtual bool isValid(const SequinID &) const = 0;
        
        static Stats analyze(Stats &, const GVariant &, const FileName &, const FileName &, const Options &, std::shared_ptr<VCFLadder>);
        
        struct TableRow
        {
            bool valid = false;
            
            Counts nr, tp, fp, fn;
            Proportion sn, pc;
            
            double fpKB;            
            Counts depth;  // Combined depth
            Counts sample; // Number of samples found
        };
        
        static TableRow getTRow(const Label &, const std::string &, const GVariant::Options &, bool keep = true);
        
        static TableRow getTotal(const Options &o) { return getTRow("NAME", "", o); }
        static TableRow getSNV(const Options &o)   { return getTRow("TYPE", "SNP", o); }
        static TableRow getIndel(const Options &o) { return getTRow("TYPE", "SNP", o, false); }
        static TableRow getHom(const Options &o)   { return getTRow("GENOTYPE", "Homozygous", o); }
        static TableRow getHet(const Options &o)   { return getTRow("GENOTYPE", "Heterozygous", o); }
        static TableRow getCode(const Options &o)  { return getTRow("GeneContext", "CodingRegion", o); }
        static TableRow getNCode(const Options &o) { return getTRow("GeneContext", "NoncodingRegion", o); }
        static TableRow getAT(const Options &o)    { return getTRow("GCcontent", "ATrich", o); }
        static TableRow getGC(const Options &o)    { return getTRow("GCcontent", "GCrich", o);  }
        static TableRow getSine(const Options &o)  { return getTRow("MobileElement", "SINE", o); }
        static TableRow getDNA(const Options &o)   { return getTRow("MobileElement", "DNA", o);  }
        static TableRow getLine(const Options &o)  { return getTRow("MobileElement", "LINE", o); }
        static TableRow getLTR(const Options &o)   { return getTRow("MobileElement", "LTR", o);  }
        static TableRow getDi(const Options &o)    { return getTRow("SimpleRepeat", "Di", o);   }
        static TableRow getTri(const Options &o)   { return getTRow("SimpleRepeat", "Tri", o);  }
        static TableRow getMono(const Options &o)  { return getTRow("SimpleRepeat", "Mono", o); }
        static TableRow getQuad(const Options &o)  { return getTRow("SimpleRepeat", "Quad", o); }
    };
}

#endif
