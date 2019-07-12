#ifndef __STATS_HPP__
#define __STATS_HPP__

#include <assert.h>
#include <algorithm>

namespace SS
{
    template <typename T1, typename T2> typename T1::value_type quant(const T1 &x, T2 q)
    {
        assert(q >= 0.0 && q <= 1.0);
        
        auto y = std::vector<typename T1::value_type>(x.begin(), x.end());
        auto const i = q * x.size();
        
        std::nth_element(y.begin(), y.begin() + i, y.end());
        return y.at(i);
    }

    template <typename T1, typename T2> typename T1::value_type hackQuant(const T1 &x, T2 q)
    {
        assert(q >= 0.0 && q <= 1.0);
        auto const i = q * x.size();
        return x.at(i);
    }

    template <typename T> typename T::value_type med(const T &x)
    {
        return quant(x, 0.5);
    }

    template <typename T> typename T::value_type min(const T &x)
    {
        return *(std::min_element(x.begin(), x.end()));
    }
    
    template <typename T> typename T::value_type max(const T &x)
    {
        return *(std::max_element(x.begin(), x.end()));
    }    
}

#endif
