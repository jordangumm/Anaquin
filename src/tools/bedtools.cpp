#include <set>
#include <iostream>
#include <fstream>
#include "data/locus.hpp"
#include "tools/tools.hpp"
#include <boost/format.hpp>
#include "tools/bedtools.hpp"
#include "parsers/parser_bed.hpp"

using namespace Anaquin;

template <typename T, typename R, template <typename, typename = std::allocator<T>> class Inputs>
static std::vector<R> inter(const Inputs<T> &x1, const Inputs<T> &x2)
{
    auto s1 = x1;
    auto s2 = x2;
    
    std::sort(s1.begin(), s1.end(), [&](const Locus &x, const Locus &y)
    {
        return (x.start < y.start) || (x.start == y.start && x.end < y.end);
    });

    std::sort(s2.begin(), s2.end(), [&](const Locus &x, const Locus &y)
    {
        return (x.start < y.start) || (x.start == y.start && x.end < y.end);
    });

    std::vector<R> r;
    
    for (auto &i : s2)
    {
        for (const auto &j : s1)
        {
            if (i.l.overlap(j.l))
            {
                auto t = i.l;
                t.inter(static_cast<Locus>(j));
                r.push_back(t);
            }
        }
    }
    
    return r;
}

template <typename T, template <typename, typename = std::allocator<T>> class Inputs>
static std::vector<std::pair<Locus, std::string>> inter2(const Inputs<T> &x1, const Inputs<T> &x2)
{
    auto s1 = x1;
    auto s2 = x2;
    
    std::sort(s1.begin(), s1.end(), [&](const Locus &x, const Locus &y)
    {
        return (x.start < y.start) || (x.start == y.start && x.end < y.end);
    });

    std::sort(s2.begin(), s2.end(), [&](const Locus &x, const Locus &y)
    {
        return (x.start < y.start) || (x.start == y.start && x.end < y.end);
    });

    std::vector<std::pair<Locus, std::string>> r;
    
    for (auto &i : s2)
    {
        for (const auto &j : s1)
        {
            if (i.l.overlap(j.l))
            {
                auto t = i.l;
                t.inter(static_cast<Locus>(j));
                r.push_back(std::pair<Locus, std::string>(t, i.name));
            }
        }
    }
    
    return r;
}

FileName BedTools::intersect2(const FileName &x, const FileName &y)
{
    std::set<ChrID> cIDs;
    std::map<ChrID, std::vector<ParserBed::Data>> r1;
    std::map<ChrID, std::vector<ParserBed::Data>> r2;
    
    ParserBed::parse(x, [&](const ParserBed::Data &i, Progress)
    {
        cIDs.insert(i.cID);
        r1[i.cID].push_back(i);
    });

    ParserBed::parse(y, [&](const ParserBed::Data &i, Progress)
    {
        cIDs.insert(i.cID);
        r2[i.cID].push_back(i);
    });

    FileName tmp = tmpFile();
    std::ofstream out(tmp);
    
    for (auto &c : cIDs)
    {
        std::sort(r1[c].begin(), r1[c].end(), [](const Locus &x, const Locus &y)
        {
            return x.start < y.start;
        });

        std::sort(r2[c].begin(), r2[c].end(), [](const Locus &x, const Locus &y)
        {
            return x.start < y.start;
        });
        
        const auto m12 = inter2<ParserBed::Data>(r1[c], r2[c]);

        if (!m12.empty())
        {
            for (auto i = 0u; i < m12.size(); i++)
            {
                out << c << "\t" << m12[i].first.start-1 << "\t" << m12[i].first.end << "\t" << m12[i].second << "\n";
            }
        }
    }
    
    out.close();
    return tmp;
}

FileName BedTools::intersect(const FileName &x, const FileName &y, Base edge)
{
    std::set<ChrID> cIDs;
    std::map<ChrID, std::vector<ParserBed::Data>> m1;
    std::map<ChrID, std::vector<ParserBed::Data>> m2;
    
    /*
     * Read the sequin regions and apply edge offset
     */
    
    ParserBed::parse(x, [&](ParserBed::Data &i, Progress)
    {
        const auto str = i.l.start + edge; // Apply edge before intersection
        const auto end = i.l.end - edge;   // Apply edge before intersection

        if (str > end)
        {
            throw std::runtime_error((boost::format("%1%:%2%-%3% is too small for edge %4%.") % i.cID % i.l.start % i.l.end % edge).str());
        }

        i.l.start = str;
        i.l.end   = end;

        cIDs.insert(i.cID);
        m1[i.cID].push_back(i);
    });

    ParserBed::parse(y, [&](const ParserBed::Data &i, Progress)
    {
        cIDs.insert(i.cID);
        m2[i.cID].push_back(i);
    });
    
    FileName tmp = tmpFile();
    std::ofstream out(tmp);

    for (auto &c : cIDs)
    {
        const auto m12 = inter<ParserBed::Data, Locus>(m1[c], m2[c]);
        
        for (auto i = 0u; i < m12.size(); i++)
        {
            // Sequin region? (e.g. GS_010)
            const auto r = std::find_if(m1[c].begin(), m1[c].end(), [&](const ParserBed::Data &x)
            {
                return x.l.overlap(m12[i]);
            });
            
            // Contig name
            const auto cont = m12.size() > 1 ? "_" + std::to_string(i+1) : "";
            
            out << c << "\t" << m12[i].start-1 << "\t" << m12[i].end << "\t" << r->name << cont << "\n";
        }
    }
    
    out.close();
    return tmp;
}
