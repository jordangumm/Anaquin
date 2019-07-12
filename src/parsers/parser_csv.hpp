#ifndef PARSER_CSV_HPP
#define PARSER_CSV_HPP

#include <vector>
#include "data/reader.hpp"

namespace Anaquin
{
    struct ParserCSV
    {
        typedef std::vector<std::string> Data;

        template <typename F> static void parse(const Reader &r, F f, const std::string &d = ",")
        {
            Progress i = 0;
            std::vector<std::string> toks;
            
            while (r.nextTokens(toks, d))
            {
                f(toks, i++);
            }
        }
    };
}

#endif
