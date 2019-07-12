#ifndef G_GERM_HPP
#define G_GERM_HPP

#include "Genomics/g_variant.hpp"

namespace Anaquin
{
    struct GGerm : public GVariant
    {
        bool isValid(const SequinID &x) const override { return isGerm(x); }

        static void report(const FileName &, const FileName &, const Options &);
    };
}

#endif
