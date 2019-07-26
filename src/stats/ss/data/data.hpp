#ifndef SS_DATA_HPP
#define SS_DATA_HPP

#include <cmath>
#include "stats/ss/data/errors.hpp"

namespace SS
{
    typedef double Real;
    
    typedef unsigned long Counts;

    typedef unsigned Level;
    typedef unsigned Factor;
    
    typedef std::size_t Index;
    
    typedef Real DF;
    typedef Real Rank;
    
    struct P
    {
        P(Real p = 0.0) : p(p)
        {
            if ((p < 0.0 || p > 1.0) && !std::isnan(p))
            {
                throw std::runtime_error("Invalid probability");
            }
        }

        P(const P &p) : p(p.p) {}

        inline Real comp() const
        {
            return P(1.0 - p);
        }

        inline void operator=(const P &p)
        {
            this->p = p.p;
        }

        // Implicity convert the probability to its underlying type
        inline operator double() const { return p; }

        Real p;
    };

    inline P operator-(Real x, const P &y)
    {
        return P(x - y.p);
    }
}

#endif
