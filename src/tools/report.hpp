#ifndef REPORT_HPP
#define REPORT_HPP

#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct Report
    {
        struct Options : public AnalyzerOptions<HumanAssembly>
        {
            Options(const AnalyzerOptions &o) : AnalyzerOptions(o) {}

            struct VOptions
            {
                // Empty Implementation
            };
            
            struct KMOptions
            {
                // Empty Implementation
            };
            
            // Generate from VCFs
            std::shared_ptr<VOptions> v;
            
            // Generate from k-mers
            std::shared_ptr<KMOptions> k;
        };

        static void rna(const Options &);
        static void meta(const Options &);
        static void genome(const Options &);
        static void somatic(const Options &);
        static void germline(const Options &);
        static void calibrate(const Options &);
    };
}

#endif
