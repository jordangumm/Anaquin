#ifndef G_SPLIT_HPP
#define G_SPLIT_HPP

#include "data/split.hpp"

namespace Anaquin
{
    struct GSplit
    {
        struct Stats : public SStats
        {
            // Median k-mer coverage for reference and variants
            std::map<SequinID, KMCoverage> VR, VV;

            // Ladder for allele frequency and synthetic
            SequinStats af, ld, gm;
        };
        
        struct Options : public SOptions
        {
            Options() : SOptions()
            {
                prod = Product::Genomics;
            }
        };

        // Build synthetic ladder
        static void buildSL(Stats &, const Options &);

        // Build allele frequency ladder
        static void buildAF(Stats &, const Options &);

        static void writeGReport(const FileName &, const SOptions &);
        static void writeQuin(const FileName &, const Stats &, const SOptions &);
        
        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
