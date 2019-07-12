#ifndef SS_DIST_HPP
#define SS_DIST_HPP

#include <math.h>
#include <ss/data/data.hpp>
#include <ss/data/errors.hpp>
#include <ss/internal/dist.hpp>

namespace SS
{
    inline double pt(double x, double n)
    {
        return Internal::pt(x, n);
    }
    
    inline double qt(double p, double ndf)
    {
        return Internal::qt(p, ndf);
    }

    inline double pnorm(double x, double mu, double sigma)
    {
        return Internal::pnorm(x, mu, sigma);
    }
}

#endif
