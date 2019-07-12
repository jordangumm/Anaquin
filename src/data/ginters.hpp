#ifndef GINTERS_HPP
#define GINTERS_HPP

#include <set>
#include <map>
#include <numeric>
#include "data/data.hpp"
#include "data/itree.hpp"
#include "data/locus.hpp"
#include "tools/errors.hpp"

namespace Anaquin
{
    class GInterval : public Matched
    {
        public:
            struct Stats
            {
                Base length = 0;
                Base nonZeros = 0;
            
                inline Proportion covered() const { return static_cast<double>(nonZeros) / length; }
            };

            GInterval() {}
            GInterval(const ChrID &cID, const Name &name, const Locus &l) : _l(l), _cID(cID), _name(name) {}

            // Return loci where no alignment
            std::set<Locus> zeros() const;
        
            Base map(const Locus &l, Base *lp = nullptr, Base *rp = nullptr);
        
            template <typename F> Stats stats(F) const
            {
                Stats stats;

                for (const auto &i : _data)
                {
                    stats.nonZeros += i.second.length();
                }
                
                stats.length = _l.length();
                assert(stats.length >= stats.nonZeros);

                return stats;
            }
        
            inline Stats stats() const
            {
                return stats([&](const ChrID &, Base, Base, Coverage)
                {
                    return true;
                });
            }

            inline void merge(const Locus &l) { _l.merge(l); }

            // Unique identifier for this interval
            inline std::string key() const { return _cID + ";" + l().key(); }
        
            inline const Locus &l() const { return _l;  }

            inline Name name() const override { return _name; }
        
        private:
        
            // The first base of the mini-region
            Base _x;
        
            // The last base of the mini-region
            Base _y;
        
            Locus _l;

            std::map<Base, Locus> _data;
        
            ChrID _cID;
        
            Name _name;
    };
    
    template <typename T = GInterval> class GIntervals
    {
        public:
        
            struct Stats : public T::Stats
            {
                // Total number of intervals
                Counts n = 0;
                
                // Number of intervals with full coverage
                Counts f = 0;
            };
        
            typedef std::map<Name, T> IntervalData;

            inline void add(const T &i)
            {
                if (_inters.count(i.key()))
                {
                    throw std::runtime_error("Duplicated key: " + i.key());
                }
                
                _inters.insert(typename std::map<Name, T>::value_type(i.key(), i));
            }

            /*
             * Merge the new interval with the first existing overlapping interval. New interval is
             * added if no overlapping found.
             */
            
            inline void merge(const T &i)
            {
                for (auto &j : _inters)
                {
                    if (j.second.l().overlap(i.l()))
                    {
                        j.second.merge(i.l());
                        return;
                    }
                }

                add(i);
            }
        
            inline void build()
            {
                std::vector<Interval_<T *>> loci;
            
                #define LOCUS_TO_TINTERVAL(x) Interval_<T *>(x.l().start, x.l().end, &x)
            
                for (auto &i : _inters)
                {
                    loci.push_back(LOCUS_TO_TINTERVAL(i.second));
                }
                
                A_CHECK(!loci.empty(), "No interval was built. Zero interval.");
            
                _tree = std::shared_ptr<IntervalTree<T *>>(new IntervalTree<T *> { loci });

                A_CHECK(_tree, "Failed to build interval treee");
            }
        
            inline T * find(const Name &x)
            {
                return _inters.count(x) ? &(_inters.at(x)) : nullptr;
            }
        
            inline const T * find(const Name &x) const
            {
                return _inters.count(x) ? &(_inters.at(x)) : nullptr;
            }

            inline T * exact(const Locus &l, std::vector<T *> *r = nullptr) const
            {
                // This could happen for chrM (no intron)
                if (!_tree)
                {
                    return nullptr;
                }
                
                auto v = _tree->findContains(l.start, l.end);

                T *t = nullptr;
    
                for (const auto &i : v)
                {
                    if (i.value->l() == l)
                    {
                        t = i.value;
                        
                        if (r)
                        {
                            r->push_back(t);
                        }
                    }
                }
            
                return t;
            }
        
            inline T * contains(const Locus &l, std::vector<T *> *r = nullptr) const
            {
                // This could happen for chrM (no intron)
                if (!_tree)
                {
                    return nullptr;
                }

                auto v = _tree->findContains(l.start, l.end);

                if (r)
                {
                    for (const auto &i : v)
                    {
                        if (i.value)
                        
                        r->push_back(i.value);
                    }
                }
            
                return v.empty() ? nullptr : v.front().value;
            }
        
            inline T * overlap(const Locus &l, std::vector<T *> *r = nullptr) const
            {
                // This could happen for chrM (no intron)
                if (!_tree)
                {
                    return nullptr;
                }

                auto v = _tree->findOverlapping(l.start, l.end);
            
                if (r)
                {
                    for (const auto &i : v)
                    {
                        r->push_back(i.value);
                    }
                }
            
                return v.empty() ? nullptr : v.front().value;
            }

            typename GIntervals::Stats stats() const
            {
                GIntervals::Stats stats;
            
                for (const auto &i : _inters)
                {
                    const auto s = i.second.stats();
                
                    stats.n++;
                    stats.length   += s.length;
                    stats.nonZeros += s.nonZeros;
                    
                    assert(s.length >= s.nonZeros);
                    
                    if (s.length == s.nonZeros)
                    {
                        stats.f++;
                    }
                }

                assert(stats.n >= stats.f);
                return stats;
            }
        
            inline const IntervalData &data() const { return _inters; }
        
            inline Base length(const Locus &l) const
            {
                std::vector<T *> o;
                overlap(l, &o);
                
                return std::accumulate(o.begin(), o.end(), 0, [&](int sums, const T *x)
                {
                    return sums + x->l().length();
                });
            }

        private:
        
            std::shared_ptr<IntervalTree<T *>> _tree;
            IntervalData _inters;
    };
}

#endif
