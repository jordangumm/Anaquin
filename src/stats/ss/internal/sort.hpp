#ifndef SS_INTERNAL3_HPP
#define SS_INTERNAL3_HPP

#include <vector>
#include <numeric>
#include <algorithm>
#include <ss/data/data.hpp>

namespace SS
{
    namespace Internal
    {
        typedef std::vector<std::size_t> Permutation;
        
        template <typename T, typename Compare> Permutation getPerm(const T& x, Compare comp)
        {
            std::vector<std::size_t> p(x.size());
            std::iota(p.begin(), p.end(), 0);
            std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) { return comp(x[i], x[j]); });
            return p;
        }
        
        template <typename T> Permutation getPerm(const T& x)
        {
            return getPerm(x, [&](typename T::value_type x, typename T::value_type y) { return x < y; });
        }

        template <typename T> std::vector<T> applyPerm(const std::vector<T>& vec, const Permutation &p)
        {
            std::vector<T> sorted(p.size());
            std::transform(p.begin(), p.end(), sorted.begin(), [&](std::size_t i){ return vec[i]; });
            return sorted;
        }
        
        template <typename T1, typename T2, typename Compare> void permSort(std::vector<T1> &x, std::vector<T2> &y, Compare comp)
        {
            auto p = getPerm(x, comp);
            
            x = applyPerm(x, p);
            y = applyPerm(y, p);
        }
        
        template <typename T1, typename T2, typename Compare> void permSort(const std::vector<T1> &x, std::vector<T2> &y, Compare comp)
        {
            auto p = getPerm(x, comp);
            y = applyPerm(y, p);
        }
        
        template <typename T1, typename T2> void permSort(std::vector<T1> &x, std::vector<T2> &y)
        {
            permSort(x, y, [&](T1 a, T1 b) { return a < b; });
        }
        
        template <typename T1, typename T2> void permSort(const std::vector<T1> &x, std::vector<T2> &y)
        {
            permSort(x, y, [&](T1 a, T1 b) { return a < b; });
        }
    }
}

#endif
