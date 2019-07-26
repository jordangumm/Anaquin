#ifndef INTERVALS_HPP
#define INTERVALS_HPP

#include <map>
#include <numeric>
#include "data/data.hpp"
#include "data/itree.hpp"
#include "data/locus.hpp"
#include "stats/stats.hpp"

namespace Anaquin
{
    class DInter : public Matched
    {
        public:
        
            typedef std::string IntervalID;
        
            struct Stats
            {
                Coverage min = std::numeric_limits<Coverage>::max();
                Coverage max = std::numeric_limits<Coverage>::min();

                // Percentile for the interval
                Coverage p25, p50, p75;
            
                // Arithmetic first moment
                double mean;
            
                // Sum of all coverage
                Coverage sums = 0;
            
                // Total number of alignments
                Counts aligns = 0;
                
                // Distribution for the coverage
                std::map<Coverage, Counts> hist;
            
                // Length of the interval
                Base length = 0;
            
                // Number of bases with non-zero coverage
                Base nonZeros = 0;
            
                // Number of bases with zero coverage
                Counts zeros = 0;
            
                std::vector<Coverage> raws;
                
                inline Proportion covered() const { return static_cast<double>(nonZeros) / length; }
            };
        
            DInter(const IntervalID &id, const Locus &l) : _id(id), _l(l)
            {
                _covs.resize(l.length());
            }

            inline void add(const Locus &l)
            {
                if ((std::size_t) l.start < _covs.size()) { _covs[l.start].starts++; }
                if ((std::size_t) l.end   < _covs.size()) { _covs[l.end].ends++;     }
                else                        { _covs.back().ends++;     }
            };

            inline Base map(const Locus &l1, Base *lp = nullptr, Base *rp = nullptr)
            {
                if (l1.start > _l.end || l1.end < _l.start)
                {
                    return 0;
                }
                
                const auto l2    = Locus(std::max(l1.start, _l.start), std::min(l1.end, _l.end));
                const auto start = l2.start - _l.start;
                const auto end   = l2.end - _l.start;

                /*
                 * For example, if the interval is (2,2038) and the locus is (2028,2042).
                 *
                 *     start = 2028 -> 2026
                 *     end   = 2038 -> 2036
                 */
            
                if (start <= end)
                {
                    _covs[start].starts++;
                    _covs[end].ends++;
                    assert(start < _covs.size() && end < _covs.size());
                }
            
                // Bases to the left of the interval fails to map
                const auto left = ((l1.start < _l.start) ? _l.start - l1.start : 0);
            
                // Bases to the right of the interval fails to map
                const auto right = ((l1.end > _l.end) ?  l1.end - _l.end : 0);
            
                if (lp) { *lp = left;  }
                if (rp) { *rp = right; }
            
                _counts++;
                
                return left + right;
            }
        
            template <typename F> Stats stats(F f) const
            {
                Stats stats;
                
                bedGraph([&](const ChrID &id, Base i, Base j, Coverage cov)
                {
                    // Should this be counted? For example, aligning to sequins?
                    if (!f(id, i, j, cov))
                    {
                        return;
                    }

                    stats.min = std::min(stats.min, cov);
                    stats.max = std::max(stats.max, cov);
                    
                    // The interval is half-open
                    const auto n = j - i;
                    
                    stats.sums      += (n * cov);
                    stats.length    += n;
                    stats.hist[cov] += n;
                    
                    for (auto i = 0; i < n; i++)
                    {
                        stats.raws.push_back(cov);
                    }
                    
                    if (!cov) { stats.zeros    += n; }
                    else      { stats.nonZeros += n; }
                });
            
                std::sort(stats.raws.begin(), stats.raws.end());
                
                stats.mean   = stats.sums / stats.length;
                stats.p25    = quant(stats.raws, 0.25);
                stats.p50    = quant(stats.raws, 0.50);
                stats.p75    = quant(stats.raws, 0.75);
                stats.aligns = count();

                return stats;
            }
        
            inline Stats stats() const
            {
                return stats([&](const ChrID, Base, Base, Coverage)
                {
                    return true;
                });
            }
        
            template <typename T> void bedGraph(T t) const
            {
                Base depth = 0;
                long lastStart = -1;
                long lastDepth = -1;
            
                for (auto j = 0; j < _l.length(); j++)
                {
                    depth += _covs[j].starts;
                
                    if (depth != lastDepth)
                    {
                        if (lastDepth != -1)
                        {
                            t(_id, lastStart, j, lastDepth);
                        }
                    
                        // Set current position as the new interval start + depth
                        lastStart = j;
                        lastDepth = depth;
                    }
                
                    depth = depth - _covs[j].ends;
                }
            
                // Print information about the last position
                if (lastDepth != -1)
                {
                    t(_id, lastStart, _l.length(), lastDepth);
                }
            }

            inline const Locus &l()       const { return _l;  }
            inline const IntervalID &id() const { return _id; }
        
            inline IntervalID name() const override { return id(); }
        
            // Return total number of alignments mapped
            inline Counts count() const { return _counts; }

        private:
        
            struct Depth
            {
                Base starts;
                Base ends;
            };
        
            IntervalID _id;
        
            // The represented interval
            Locus _l;
        
            // Number of alignments mapped to the interval
            Counts _counts = 0;

            // For each base in the interval (relative to the beginning of the interval)
            std::vector<Depth> _covs;
    };
    
    template <typename T = DInter> class DIntervals
    {
        public:
        
            typedef std::map<typename T::IntervalID, T> IntervalData;

            inline void add(const T &i)
            {
                _inters.insert(typename std::map<typename T::IntervalID, T>::value_type(i.id(), i));
            }

            inline void build()
            {
                std::vector<Interval_<T *>> loci;
            
                #define LOCUS_TO_TINTERVAL(x) Interval_<T *>(x.l().start, x.l().end, &x)
            
                for (auto &i : _inters)
                {
                    loci.push_back(LOCUS_TO_TINTERVAL(i.second));
                }
                
                if (loci.empty())
                {
                    throw std::runtime_error("No interval was built. loci.empty().");
                }
            
                _tree = std::shared_ptr<IntervalTree<T *>>(new IntervalTree<T *> { loci });
            }
        
            inline T * find(const typename T::IntervalID &id)
            {
                return _inters.count(id) ? &(_inters.at(id)) : nullptr;
            }
        
            inline const T * find(const typename T::IntervalID &id) const
            {
                return _inters.count(id) ? &(_inters.at(id)) : nullptr;
            }

            inline T * exact(const Locus &l, std::vector<T *> *r = nullptr) const
            {
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

            template <typename F> void bedGraph(F f) const
            {
                for (const auto &i : _inters)
                {
                    i.second.bedGraph([&](const ChrID &id, Base i, Base j, Coverage cov)
                    {
                        f(id, i, j, cov);
                    });
                }
            }

            typename T::Stats stats() const
            {
                DInter::Stats stats;
            
                for (const auto &i : _inters)
                {
                    const auto s = i.second.stats();
                
                    stats.sums     += s.sums;
                    stats.aligns   += s.aligns;
                    stats.length   += s.length;
                    stats.nonZeros += s.nonZeros;
                    stats.zeros    += s.zeros;
                    stats.min       = std::min(stats.min, s.min);
                    stats.max       = std::max(stats.max, s.max);
                
                    for (const auto &j : s.hist)
                    {
                        stats.hist[j.first] += j.second;
                    }
                }
            
                stats.mean = stats.sums / stats.length;
            
                return stats;
            }
        
            inline const IntervalData &data() const { return _inters; }
        
            // Number of intervals
            inline Counts size() const { return _inters.size(); }
        
            inline Base length() const
            {
                return std::accumulate(_inters.begin(), _inters.end(), 0,
                        [&](int sums, const std::pair<DInter::IntervalID, DInter> & p)
                {
                    return sums + p.second.l().length();
                });
            }
        
        private:
        
            std::shared_ptr<IntervalTree<T *>> _tree;
        
            IntervalData _inters;
    };

    typedef std::map<ChrID, DIntervals<>> Chr2DInters;

    template <typename T, typename D = DInter> D * contains(const T &x, const ChrID &cID, const Locus &l)
    {
        if (!x.count(cID))
        {
            return nullptr;
        }
  
        return x.at(cID).contains(l);
    }
    
    struct ID2Intervals : std::map<DInter::IntervalID, DIntervals<>>
    {
        inline void add(const DInter::IntervalID &id, const DIntervals<> &i)
        {
            (*this)[id] = i;
        }
        
        template <typename F> void bedGraph(F f) const
        {
            for (const auto &i : (*this))
            {
                i.second.bedGraph([&](const ChrID &id, Base i, Base j, Coverage cov)
                {
                    f(id, i, j, cov);
                });
            }
        }

        inline DInter::Stats stats() const
        {
            DInter::Stats stats;
            
            for (const auto &i : *this)
            {
                const auto s = i.second.stats();
                
                stats.sums     += s.sums;
                stats.length   += s.length;
                stats.nonZeros += s.nonZeros;
                stats.zeros    += s.zeros;
                stats.aligns   += s.aligns;
                stats.min       = std::min(stats.min, s.min);
                stats.max       = std::max(stats.max, s.max);
                
                for (const auto &j : s.hist)
                {
                    stats.hist[j.first] += j.second;
                }
            }
            
            auto percent = [&](Counts n)
            {
                Counts i = 0;
                
                for (const auto &p : stats.hist)
                {
                    i += p.second;
                    
                    // Have we reached our percentile?
                    if (i >= n)
                    {
                        return p.first;
                    }
                }
                
                return stats.max;
            };
            
            stats.mean = stats.sums / stats.length;
            stats.p25  = percent(0.25 * stats.length);
            stats.p50  = percent(0.50 * stats.length);
            stats.p75  = percent(0.75 * stats.length);

            return stats;
        }
    };
}

#endif
