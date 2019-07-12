#ifndef R_SPLIT_HPP
#define R_SPLIT_HPP

#include "data/split.hpp"
#include "stats/linear.hpp"

namespace Anaquin
{
    struct RSplit
    {
        struct Stats
        {
            // Before calibration
            SStats S1;
            
            // After percentage calibration for "sequins"
            SStats S2;
            
            // Calibration statistics for "sequins"
            SCStats C1;
            
            // Isoform and gene-level expression
            SequinStats l1, l2;
        };
        
        struct Options : public SOptions
        {
            Options()
            {
                prod = Product::RNA;
            }
            
            Mixture mix;
        };
        
        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
