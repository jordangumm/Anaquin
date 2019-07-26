#ifndef G_SYNTHETIC_HPP
#define G_SYNTHETIC_HPP

#include "data/data.hpp"
#include "Kallisto.hpp"

namespace Anaquin
{
    struct LadderInternal
    {
        struct Results
        {
            // Average scaling factor (different for each ladder standard)
            Proportion meanS;
            
            FileName r1, r2;
            
            // Coverage for each ladder standard
            std::map<StandardID, KMCoverage> cov;

            // Counts for the target
            Counts targetC = 0;
            
            // Calibrated counts
            Counts after = 0;
        };
        
        static Results calibrate(const KStats &,
                                 const Path &,
                                 const KOptions &,
                                 Proportion,
                                 const WriterOptions &,
                                 const Label &);
    };
}

#endif
