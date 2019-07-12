#include <data/ginters.hpp>

using namespace Anaquin;

std::set<Locus> GInterval::zeros() const
{
    std::set<Locus> r;
    
    Base i = 1;
    
    for (auto &j : _data)
    {
        if (i < j.second.start)
        {
            r.insert(Locus(i, j.second.start-1));
        }
        
        i = j.second.end + 1;
    }

    return r;
}

Base GInterval::map(const Locus &l, Base *lp, Base *rp)
{
    bool added = true;
    bool p1 = true;
    bool p2 = false;
    
    // Pointing to the matching
    Locus *m = nullptr;
    
    Base left  = 0;
    Base right = 0;
    
    Base origKey = 0;
    
    auto fixKey = [&]()
    {
        assert(_data.count(origKey));
        assert(m && origKey && m->start && m->end);
        
        auto iter = _data.find(origKey);
        assert(iter != _data.end());
        
        // Remove the old key
        _data.erase(iter);
        
        auto l = *m;
        
        // Add the new key
        _data[l.end] = l;
    };
    
    if (_x == _y && _x == 0)
    {
        added = true;
    }
    else if (l.end < _x || l.start > _y)
    {
        added = true;
    }
    else
    {
        // Can we quickly find a match by log(n)?
        auto iter = _data.lower_bound(l.start);
        
        auto it = iter;

        if (!it->second.overlap(l))
        {
            added = true;
        }
        else
        {
            assert(it != _data.end());
            assert(it->second.overlap(l));
            
            for (; it != _data.cend();)
            {
                auto &j = it->second;
                
                if (!j.overlap(l))
                {
                    break;
                }
                
                if (p1)
                {
                    if (j.overlap(l))
                    {
                        origKey = it->first;
                        added = false;
                        
                        left  = ((l.start < j.start) ? j.start -  l.start : 0);
                        right = ((l.end   > j.end)   ?  l.end  - j.end   : 0);
                        
                        j.start = std::min(j.start, l.start);
                        j.end   = std::max(j.end,   l.end);
                        j.start = std::max(j.start, _l.start);
                        j.end   = std::min(j.end,   _l.end);

                        // So that we can access it in the later iterations
                        m = &j;
                        
                        if (lp) { *lp = left;  }
                        if (rp) { *rp = right; }
                        
                        p1 = false;
                        p2 = true;
                    }
                }
                else if (p2)
                {
                    if (!j.overlap(l))
                    {
                        _x = _data.begin()->second.start;
                        _y = _data.end()->second.end;
                        
                        fixKey();

                        return (left + right);
                    }
                    
                    m->start = std::min(m->start, j.start);
                    m->end   = std::max(m->end,   j.end);
                    m->start = std::max(m->start, _l.start);
                    m->end   = std::min(m->end,   _l.end);
                    
                    // Remove the overlapping entry
                    _data.erase(it++);
                    
                    continue;
                }
                else
                {
                    throw std::runtime_error("Invalid phase");
                }
                
                ++it;
            }
        }
    }
    
    if (added)
    {
        auto start = l.start;
        auto end = l.end;
        
        start = std::max(start, _l.start);
        end   = std::min(end,   _l.end);
        
        _data[end] = Locus(start, end);
        
        left  = ((l.start < _l.start) ? _l.start -  l.start : 0);
        right = ((l.end   > _l.end)   ?  l.end  - _l.end   : 0);
    }
    else
    {
        fixKey();
    }
    
    if (lp) { *lp = left;  }
    if (rp) { *rp = right; }
    
    _x = _data.begin()->second.start;
    _y = _data.rbegin()->second.end;
    
    return left + right;
}
