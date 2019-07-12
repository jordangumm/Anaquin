#ifndef BEDTOOLS_HPP
#define BEDTOOLS_HPP

#include "data/data.hpp"

namespace Anaquin
{
    struct BedTools
    {
        static FileName intersect(const FileName &, const FileName &, Base edge);
        static FileName intersect2(const FileName &, const FileName &);
    };
}

#endif
