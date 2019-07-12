#ifndef B_DATA_HPP
#define B_DATA_HPP

#include "tools/tools.hpp"
#include "data/ginters.hpp"
#include "data/dinters.hpp"
#include <boost/format.hpp>
#include "parsers/parser_bed.hpp"

namespace Anaquin
{
    typedef ParserBed::Data Data;

    struct BedChrData
    {
        std::map<Locus, Data> l2d;
    };
    
    struct BedData : public std::map<ChrID, BedChrData>
    {
        inline std::set<Name> names() const
        {
            std::set<Name> x;
            
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.l2d)
                {
                    x.insert(j.second.name);
                }
            }
            
            return x;
        }

        inline const ParserBed::Data * find(const Name &x)
        {
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.l2d)
                {
                    if (j.second.name == x)
                    {
                        return &j.second;
                    }
                }
            }
            
            return nullptr;
        }

        inline Base len(const ChrID &x, const Name &y = "") const
        {
            Base b = 0;
            
            for (const auto &i : at(x).l2d)
            {
                if (y.empty() || i.second.name == y)
                {
                    b += i.second.l.length();
                }
            }
            
            return b;
        }

        inline Base len() const
        {
            return countMap(*this, [&](const ChrID &x, const BedChrData &)
            {
                return len(x);
            });
        }
        
        inline Counts count() const
        {
            return countMap(*this, [&](const ChrID &, const BedChrData &x)
            {
                return x.l2d.size();
            });
        }
        
        inline DIntervals<> inters(const ChrID &x) const
        {
            DIntervals<> r;
            
            for (const auto &i : at(x).l2d)
            {
                r.add(DInter(i.second.name, i.second.l));
            }
            
            r.build();
            return r;
        }

        inline std::map<ChrID, DIntervals<>> inters() const
        {
            std::map<ChrID, DIntervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = inters(i.first);
            }
            
            return r;
        }
        
        inline const ParserBed::Data * overlap(const ChrID &x, const Locus &l)
        {
            for (const auto &i : *this)
            {
                if (i.first == x)
                {
                    for (const auto &j : i.second.l2d)
                    {
                        if (j.first.overlap(l))
                        {
                            return &j.second;
                        }
                    }
                }
            }
            
            return nullptr;
        }

        inline GIntervals<> ginters(const ChrID &x) const
        {
            GIntervals<> r;
            
            for (const auto &i : at(x).l2d)
            {
                r.add(GInterval(x, i.second.name, i.second.l));
            }
            
            r.build();
            return r;
        }
        
        inline std::map<ChrID, GIntervals<>> ginters() const
        {
            std::map<ChrID, GIntervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = ginters(i.first);
            }
            
            return r;
        }
        
        // Source for regions (if available)
        FileName src;
    };
    
    struct RegionOptions
    {
        RegionOptions(Base edge = 0) : edge(edge) {}
        
        Base edge = 0;
        
        // Only those sequins?
        std::set<SequinID> only;
        
        // Any sequin containing the substrings?
        std::set<std::string> conts;
    };
    
    template <typename F> BedData readRegions(const Reader &r, F f, RegionOptions o = RegionOptions())
    {
        BedData c2d;
        
        ParserBed::parse(r, [&](ParserBed::Data &x, Progress i)
        {
            if (x.l.length() < 2 * o.edge)
            {
                throw std::runtime_error((boost::format("%1% %2%-%3% is too narrow. Please check the edge width.")
                                                % x.cID
                                                % x.l.start
                                                % x.l.end).str());
            }
            
            /*
             * Filter out unwanted regions if we have the information (fourth column in BED). Otherwise, don't do it.
             */
            
            else if (!o.only.empty() && !o.only.count(x.name) && !std::any_of(o.conts.begin(), o.conts.end(), [&](const SequinID &i)
            {
                return isSubstr(x.name, i);
            }))
            {
                return;
            }
            
            x.l.end   -= o.edge;
            x.l.start += o.edge;

            c2d[x.cID].l2d[x.l] = x;            
            f(x, i);
        });

        c2d.src = r.src();
        assert(!c2d.src.empty());
        
        return c2d;
    }
    
    inline BedData readRegions(const Reader &r)
    {
        return readRegions(r, [](const ParserBed::Data &, Progress) {});
    }
}

#endif
