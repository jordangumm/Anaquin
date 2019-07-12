#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <functional>
#include <klib/khash.h>
#include "stats/analyzer.hpp"
#include "writers/sam_writer.hpp"

namespace Anaquin
{
    struct Alignment;

    class Random
    {
        public:
        
            Random(double prob) : _prob(prob)
            {
                assert(prob >= 0.0);
                _seed = rand();
            }

            inline bool select(const std::string &hash) const
            {
                const uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(hash.c_str()) ^ _seed);
                return ((double)(k&0xffffff) / 0x1000000 >= _prob);
            }
        
        private:
        
            // Random seed
            int _seed;
        
            // The probability of selection
            Probability _prob;
    };

    struct Sampler
    {
        struct SGReads
        {
            Reads syn = 0, gen = 0;
            
            inline Proportion dilut() const
            {
                return static_cast<Proportion>(syn) / (syn + gen);
            }
        };

        struct Stats
        {
            SGReads before, after;
        };
        
        template <typename O> static Stats sample(const FileName &file,
                                                  Proportion p,
                                                  const O &o,
                                                  std::function<bool (const ChrID &)> isSyn)
        {
            Sampler::Stats stats;
            
            A_ASSERT(p > 0.0 && p <= 1.0);
            Random r(1.0 - p);
            
            SAMWriter w(true, true);
            w.open("");
            
            ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &info)
            {
                if (info.p.i && !(info.p.i % 1000000))
                {
                    o.logInfo(std::to_string(info.p.i));
                }
                
                const auto shouldWrite = !x.mapped || !isSyn(x.cID) || !x.isAligned;
                
                if (x.isPrimary && x.isAligned)
                {
                    if (isSyn(x.cID))
                    {
                        stats.before.syn++;
                    }
                    else
                    {
                        stats.before.gen++;
                    }
                }
                
                // This is the key, randomly write the reads with certain probability
                if (shouldWrite || r.select(x.name))
                {
                    if (x.isPrimary && x.isAligned && isSyn(x.cID))
                    {
                        stats.after.syn++;
                        o.logInfo("Sampled " + x.name);
                    }
                    
                    /*
                     * TopHat2 might give an empty QNAME, which violates the SAM/BAM format. It's fine to
                     * give '*' to QNAME, but not an empty string....
                     */
                    
                    if (!x.name.empty())
                    {
                        // Print SAM line
                        w.write(x);
                    }
                }
            }, true);
            
            A_ASSERT(stats.before.syn >= stats.after.syn);
            stats.after.gen = stats.before.gen;
            
            w.close();
            
            return stats;
        }
    };
}

#endif
