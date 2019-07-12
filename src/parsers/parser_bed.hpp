#ifndef PARSER_BED_HPP
#define PARSER_BED_HPP

#include "data/locus.hpp"
#include "data/reader.hpp"
#include <boost/algorithm/string.hpp>

namespace Anaquin
{
    struct ParserBed
    {
        struct Data
        {
            operator const Locus &() const { return l; }
            operator const std::string &() const { return name; }
            
            ChrID cID;
            
            Locus l;

            // 4-th column in BED
            std::string name;
        };

        template <typename F> static void parse(const Reader &r, F f)
        {
            Data d;
            Progress i = 0;
            
            std::vector<std::string> sizes, starts, tokens;
            
            while (r.nextTokens(tokens, "\t"))
            {
                // Empty line?
                if (tokens.size() == 1)
                {
                    return;
                }
                
                // Name of the chromosome
                d.cID = CHROM(tokens[0]);

                /*
                 * https://genome.ucsc.edu/FAQ/FAQformat.html#format1
                 *
                 * The starting position requires to increment because BED position is a 0-based.
                 * The end position requies no increment because it is "not included in the display".
                 */
                
                d.l = Locus(stod(tokens[1])+1, stod(tokens[2]));
                
                if (tokens.size() >= 4)
                {
                    // Name of the BED line (eg: gene)
                    d.name = tokens[3];
                }
                
                f(d, i++);
            }
        }
    };
}

#endif
