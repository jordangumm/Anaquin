#ifndef HIST_HPP
#define HIST_HPP

#include <set>
#include <map>
#include <string>
#include <assert.h>
#include "data/data.hpp"

namespace Anaquin
{
    typedef std::map<std::string, Counts> Hist;
    
    typedef std::map<SequinID, Counts> SequinHist;
    
    template <typename T> Hist createHist(const T& t)
    {
        Hist hist;
        
        for (const auto &i : t)
        {
            hist[i.first] = 0;
        }
        
        assert(!hist.empty());
        return hist;
    }
    
    template <typename T> Hist createHist(const std::set<T> & t)
    {
        Hist hist;
        
        for (const auto &i : t)
        {
            hist[i] = 0;
        }
        
        assert(!hist.empty());
        return hist;
    }
}

#endif
