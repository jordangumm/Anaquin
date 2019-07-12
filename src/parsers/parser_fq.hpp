#ifndef PARSER_FQ_HPP
#define PARSER_FQ_HPP

#include "tools/tools.hpp"
#include "data/reader.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserFQ
    {
        struct Data
        {
            Sequence seq1,  seq2;
            ReadName name1, name2;
            Sequence opt1,  opt2;
            Sequence qual1, qual2;
        };

        template <typename F> static void parse(const Reader &r1, const Reader &r2, F f)
        {
            Data x;
            std::string l1, l2;
            
            while (r1.nextLine(l1) && r2.nextLine(l2))
            {
                x.name1 = remove(l1, "@");
                x.name2 = remove(l2, "@");

                r1.nextLine(l1);
                r2.nextLine(l2);
                x.seq1 = l1;
                x.seq2 = l2;

                r1.nextLine(l1);
                r2.nextLine(l2);
                x.opt1 = l1;
                x.opt2 = l2;

                r1.nextLine(l1);
                r2.nextLine(l2);
                x.qual1 = l1;
                x.qual2 = l2;

                f(x);
            }
        }
    };
}

#endif
