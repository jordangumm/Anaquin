#ifndef G_SYNTHETIC_HPP
#define G_SYNTHETIC_HPP

#include "Kallisto.hpp"

namespace Anaquin
{
    struct GSynthetic
    {
        struct Results
        {
            FileName r1, r2;
            
            // Coverage for each synthetic standard
            std::map<StandardID, KMCoverage> cov;

            // Calibration target
            SequinID target;
            
            // Counts for the target
            Counts targetC;
        };
        
        static Results calibrate(const KStats &,
                                 const Path &,
                                 const KOptions &,
                                 const WriterOptions &,
                                 const Label &);
    };
}

#endif
