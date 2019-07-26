#ifndef BED_TOOLS_HPP
#define BED_TOOLS_HPP

#include "data/data.hpp"

namespace Anaquin
{
    struct BedTools
    {
        static FileName intersect(const FileName &, const FileName &, Base);
        static FileName intersect2(const FileName &, const FileName &);
    };
}

#endif
