#ifndef LADDER_HPP
#define LADDER_HPP

#include <map>
#include <algorithm>
#include "data/data.hpp"
#include "tools/tools.hpp"

namespace Anaquin
{
    struct Ladder
    {
        inline void add(const SequinID &x, Mixture m, Concent c)
        {
            seqs.insert(x);
            
            switch (m)
            {
                case Mix_1: { m1[x] = c; break; }
                case Mix_2: { m2[x] = c; break; }
            }
        }
        
        inline Counts count(Concent i, Mixture m)
        {
            const auto &p = m == Mix_1 ? m1 : m2;
            
            return std::count_if(p.begin(), p.end(), [&](const std::pair<SequinID, Concent> &x)
            {
                return x.second == i;
            });
        }
        
        inline Counts count() const { return seqs.size(); }

        inline bool contains(const SequinID &x, Mixture m = Mix_1)
        {
            return m == Mix_1 ? m1.count(x) : m2.count(x);
        }
        
        inline Concent input(const SequinID &x, Mixture m = Mix_1)
        {
            return m == Mix_1 ? m1.at(x) : m2.at(x);
        }

        inline bool hasMix2() const { return !m2.empty(); }
        
        FileName src;
        
        std::set<SequinID> seqs;
        
        // Ladder for mixture 1
        std::map<SequinID, Concent> m1;
        
        // Ladder for mixture 2
        std::map<SequinID, Concent> m2;
    };
}

#endif
