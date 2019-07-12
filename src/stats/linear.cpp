#include "stats/linear.hpp"
#include "tools/errors.hpp"
#include <ss/regression/linear.hpp>

using namespace Anaquin;

Limit SequinStats::limitQuant() const
{
    Limit limit;
    
    for (const auto &i : *this)
    {
        if (limit.id.empty() || limit.abund > i.second.x)
        {
            limit.id = i.first;
            limit.abund = i.second.x;
        }
    }

    return limit;
}

SequinStats::Data SequinStats::data(bool shouldLog, bool ignoreZero) const
{
    Data r;
    
    auto f = [&](double x)
    {
        return shouldLog ? (log2(x ? x : 1)) : x;
    };
    
    for (const auto &p : *this)
    {
        if (!std::isnan(p.second.x) && !std::isnan(p.second.y))
        {
            if (!ignoreZero || p.second.y != 0.0)
            {
                const auto x = f(p.second.x);
                const auto y = f(p.second.y);
                
                r.x.push_back(x);
                r.y.push_back(y);
                r.ids.push_back(p.first);
                
                r.id2x[p.first] = x;
                r.id2y[p.first] = y;                
            }
        }
    }
    
    return r;
}

LinearModel SequinStats::linear(bool shouldLog, bool ignoreZero) const
{
    const auto d = data(shouldLog, ignoreZero);
    
    LinearModel lm;
    
    try
    {
        A_CHECK(!d.x.empty() && !d.y.empty(), "Failed to perform linear regression. Empty inputs.");
        A_CHECK(std::adjacent_find(d.x.begin(), d.x.end(), std::not_equal_to<double>()) != d.x.end(),
                        "Failed to perform linear regression. Flat mixture?");

        const auto m = SS::linearModel(d.y, d.x);

        lm.F  = m.f;
        lm.p  = m.p;
        lm.r  = SS::pearson(d.x, d.y);
        lm.c  = m.coeffs[0].est;
        lm.m  = m.coeffs[1].est;
        lm.R2 = m.r2;
    }
    catch(...)
    {
        lm.F  = NAN;
        lm.p  = NAN;
        lm.r  = NAN;
        lm.c  = NAN;
        lm.m  = NAN;
        lm.R2 = NAN;
    }

    return lm;
}
