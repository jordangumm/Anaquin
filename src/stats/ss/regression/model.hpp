#ifndef SS_MODELT_HPP
#define SS_MODELT_HPP

#include "stats/ss/data/data.hpp"

namespace SS
{
    struct Coefficient
    {
        Real est;

        // Standard error of the coefficient
        Real se;

        // Test statistic whether the coefficient is significant from zero
        Real test;
        
        // P-value for the test statistic
        Real p;

        // Upper 95% confidence interval
        Real upper;
        
        // Lower 95% confidence interval
        Real lower;
        
        DF df;
    };
    
    struct Variation
    {
        // Sum-of-squares of the variation
        Real ss;
        
        // Mean square of the variation
        Real ms;

        // Degrees of freedom
        DF df;
    };
}

#endif
