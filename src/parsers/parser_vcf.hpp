#ifndef PARSER_VCF_HPP
#define PARSER_VCF_HPP

#include <map>
#include "data/reader.hpp"
#include "data/variant.hpp"

namespace Anaquin
{
    struct ParserVCF
    {
        typedef std::function<void (Variant &)> Functor;
        static void parse(const Reader &, Functor);
        
        /*
         * Parse for all contigs in VCF
         */
        
        static std::map<ChrID, Base> contigs(const Reader &);
    };
}

#endif
