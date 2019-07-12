#ifndef SS_INTERNAL_STATS_HPP
#define SS_INTERNAL_STATS_HPP

#include <numeric>
#include <algorithm>
#include <ss/data/data.hpp>

namespace SS
{
    namespace Internal
    {
        template <typename T> typename T::value_type sum(const T &x)
        {
            return std::accumulate(x.begin(), x.end(), 0.0);
        }
        
        template <typename T> Counts count(const T &x)
        {
            return static_cast<Counts>(std::distance(x.begin(), x.end()));
        }
        
        template <typename T> typename T::value_type cov(const T &x, const T &y)
        {
            const auto sum_x = std::accumulate(std::begin(x), std::end(x), 0.0);
            const auto sum_y = std::accumulate(std::begin(y), std::end(y), 0.0);
            
            const auto mx = sum_x / x.size();
            const auto my = sum_y / y.size();
            
            auto accum = 0.0;
            
            auto i = x.begin();
            auto j = y.begin();
            
            for (; i != x.end(); i++, j++)
            {
                accum += (*i - mx) * (*j - my);
            }
            
            return accum / (x.size() - 1.0);
        }
    }
}

#endif
