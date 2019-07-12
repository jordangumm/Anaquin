#ifndef CALIBRATE_HPP
#define CALIBRATE_HPP

#include "Kallisto.hpp"
#include "data/data.hpp"
#include "tools/random.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    class Calibrator
    {
        public:
            struct Result
            {
                // Number of reads selected
                Counts n = 0;
                
                FileName o1, o2;
            };

            virtual Result calibrate(const KStats &, Selection &, const WriterOptions &) = 0;
        
            /*
             * Calibrate given proability of selecting a pair of reads
             */
        
            virtual Result calibrate(const KStats &, Probability, const WriterOptions &) = 0;

            // Calibration for BAM format
            static std::shared_ptr<Calibrator> createBAM(const FileName &, const FileName &);
        
            // Calibration for FASTQ format
            static std::shared_ptr<Calibrator> createFQ(const FileName &, const FileName &,
                                                        const FileName &, const FileName &);
    };
}

#endif
