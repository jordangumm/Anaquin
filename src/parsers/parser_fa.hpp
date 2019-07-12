#ifndef PARSER_FA_HPP
#define PARSER_FA_HPP

#include "data/reader.hpp"
#include "stats/analyzer.hpp"
#include <boost/algorithm/string.hpp>

namespace Anaquin
{
    struct ParserFA
    {
        struct Data
        {
            ChrID id;
            Sequence seq;
        };

        typedef std::function<void(const Data &)> Callback;

        static void parse(const Reader &r, Callback f, const ChrID &chrID = "")
        {
            Data l;
            std::string s;
            Progress p = 0;
            
            std::stringstream ss;
            #define CALL_BACK() if (p) { l.seq = ss.str(); f(l); ss.str(""); }

            while (r.nextLine(s))
            {
                if (s[0] != '>')
                {
                    if (chrID.empty() || l.id == chrID)
                    {
                        boost::trim(s);
                        ss << s;
                    }
                }
                else
                {
                    CALL_BACK();
                    
                    l.id = s.substr(1, s.size() - 1);
                }
                
                p++;
            }
            
            CALL_BACK();
        }
        
        /*
         * Calculate size of each chromosome in FASTA
         */
        
        static void size(const Reader &r, std::map<SequinID, Base> &s)
        {
            ParserFA::parse(r, [&](const ParserFA::Data &x)
            {
                s[x.id] = x.seq.size();
            });
        }
    };
}

#endif
