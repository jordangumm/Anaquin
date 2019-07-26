#ifndef M_SPLIT_HPP
#define M_SPLIT_HPP

#include "data/split.hpp"
#include "stats/linear.hpp"

namespace Anaquin
{
    struct MSplit
    {
        struct Stats
        {
            // Before calibration
            SStats S1;
            
            // After percentage calibration for "sequins"
            SStats S2;
            
            // After synthetic calibration for "ladders"
            SStats S3;
            
            // After both synthetic and percentage calibration for "ladders"
            SStats S4;
            
            // Calibration statistics for stage one "sequins"
            SCStats C1;
            
            // Calibration statistics for stage one "ladders"
            SCStats C2;
            
            // Calibration statistics for stage two "ladders"
            LadderInternal::Results S5;
            
            SequinStats l1, l2;
        };

        struct Options : public SOptions
        {
            Options()
            {
                prod = Product::Meta;
            }
            
            Mixture mix;
        };
        
        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
