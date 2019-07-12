#ifndef PARSER_BAMBED_HPP
#define PARSER_BAMBED_HPP

#include "data/dinters.hpp"
#include "parsers/parser_bam.hpp"

namespace Anaquin
{
    struct ParserBAMBED
    {
        struct Stats
        {
            Chr2DInters inters;
            
            Counts nNA  = 0;
            Counts nMap = 0;
        };

        enum class Response
        {
            OK,
            SKIP_MATCH,
            SKIP_EVERYTHING
        };
        
        template <typename F> static ParserBAMBED::Stats parse(const FileName &file,
                                                               const Chr2DInters &r,
                                                               F f)
        {
            ParserBAMBED::Stats stats;

            // For each chromosome...
            for (const auto &i : r)
            {
                DIntervals<> x;
                
                for (const auto &inter : i.second.data())
                {
                    const auto &l = inter.second.l();
                    x.add(DInter(l.key(), l));
                }

                stats.inters[i.first] = x;
                stats.inters[i.first].build();
            }

            ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &info)
            {
                DInter *matched = nullptr;
                
                if (x.mapped && stats.inters.count(x.cID))
                {
                    matched = stats.inters[x.cID].overlap(x.l);
                }
                
                const auto r = f(x, info, matched);
                
                if (r != Response::SKIP_EVERYTHING && !x.isDuplicate)
                {
                    if (r != Response::SKIP_MATCH && matched)
                    {
                        matched->map(x.l);

                        if (x.mapped) { stats.nMap++; }
                        else          { stats.nNA++;  }
                    }
                    else
                    {
                        stats.nNA++;
                    }
                }
            });
            
            return stats;
        }
    };
}

#endif
