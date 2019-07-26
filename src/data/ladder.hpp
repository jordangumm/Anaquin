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
            mix(m)[x] = c;
        }
        
        inline std::map<SequinID, Concent> &mix(Mixture m)
        {
            switch (m)
            {
                case Mix_1: { return m1; }
                case Mix_2: { return m2; }
                case Mix_3: { return m3; }
            }
        }
        
        inline Counts count(Concent i, Mixture m)
        {
            return std::count_if(mix(m).begin(), mix(m).end(), [&](const std::pair<SequinID, Concent> &x)
            {
                return x.second == i;
            });
        }
        
        inline Counts count() const { return seqs.size(); }

        inline bool contains(const SequinID &x, Mixture m = Mix_1)
        {
            return mix(m).count(x);
        }
        
        inline Concent input(const SequinID &x, Mixture m = Mix_1)
        {
            return mix(m).at(x);
        }

        inline bool hasMix2() const { return !m2.empty(); }
        inline bool hasMix3() const { return !m3.empty(); }

        FileName src;
        
        std::set<SequinID> seqs;
        std::map<SequinID, Concent> m1, m2, m3;
    };
}

#endif
