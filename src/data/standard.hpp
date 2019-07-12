#ifndef STANDARD_HPP
#define STANDARD_HPP

#include "data/reference.hpp"

namespace Anaquin
{
    class Standard
    {
        public:

            static Standard& instance(bool reload = false)
            {
                static Standard s;
                
                // Reload the default resources
                if (reload)
                {
                    s = Standard();
                }
                
                return s;
            }

            // Add sequin regions in BED format
            static BedData readBED(const Reader &, const RegionOptions &);

           /*
            * ---------------- Metagenomcis analysis ----------------
            */
        
            Ladder readMMix(const Reader &);

            MetaRef meta;

           /*
            * ---------------- RNA analysis ----------------
            */

            Ladder readRLen(const Reader &);
            Ladder readRMix(const Reader &);
            Ladder readRGMix(const Reader &);

            RnaRef rna;

            /*
             * ---------------- Genomics analysis ----------------
             */

            // Both germline and somatic variants
            static VCFLadder addVCF(const Reader &, std::shared_ptr<BedData>);

            // Filter germline variants
            static VCFLadder addGVCF(const Reader &, std::shared_ptr<BedData>);
        
            // Filter somatic variants
            static VCFLadder addSVCF(const Reader &, std::shared_ptr<BedData>);
        
            Ladder addAF(const Reader &);
        
            GenomicsRef gen;

        private:
            Standard() {}
            Standard(Standard const&) = delete;
    };
}

#endif
