#ifndef SS_INTERNAL_DIST_HPP
#define SS_INTERNAL_DIST_HPP

#include <math.h>
#include "stats/ss/data/data.hpp"
#include "stats/ss/data/errors.hpp"

#ifndef HAVE_RMATH
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#endif

namespace SS
{
    namespace Internal
    {
        inline double pf(double q, double df1, double df2)
        {
#ifdef HAVE_RMATH
            return RMath::pt(x, n);
#else
            auto f = boost::math::fisher_f(df1, df2);
            return boost::math::cdf(f, q);
#endif
        }
        
        inline double pt(double x, double df)
        {
#ifdef HAVE_RMATH
            return RMath::pt(x, n);
#else
            auto d = boost::math::students_t(df);
            return boost::math::cdf(d, x);
#endif
        }
        
        inline double qt(double p, double df)
        {
#ifdef HAVE_RMATH
            return RMath::qt(x, n);
#else
            auto d = boost::math::students_t(df);
            return boost::math::quantile(d, p);
#endif
        }
        
        inline double pnorm(double x, double mu, double sigma)
        {
#ifdef HAVE_RMATH
            return RMath::pnorm(x, mu, sigma);
#else
            auto d = boost::math::normal(mu, sigma);
            return boost::math::cdf(d, x);
#endif
        }
    }
}

#endif
