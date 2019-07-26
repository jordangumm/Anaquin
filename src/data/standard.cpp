#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "RNA/RNA.hpp"
#include "data/reader.hpp"
#include "tools/errors.hpp"
#include "data/standard.hpp"
#include "Genomics/Genomics.hpp"
#include "parsers/parser_csv.hpp"
#include "parsers/parser_bed.hpp"

using namespace Anaquin;

enum MixtureFormat
{
    X_M,
    M_X_M,
    X_X_X_M,
    X_M_X_M,
};

static unsigned countColumns(const Reader &r)
{
    std::size_t n = 0;

    ParserCSV::parse(Reader(r), [&](const ParserCSV::Data &d, Progress)
    {
        n = std::max(n, d.size());
    }, "\t");

    return static_cast<unsigned>(n);
}

BedData Standard::readBED(const Reader &r, const RegionOptions &o)
{
    return readRegions(Reader(r), [&](const ParserBed::Data &, Progress) {}, o);
}

template <typename Reference> Ladder readAF(const Reader &r, Reference &)
{
    Ladder x = Ladder();
    
    ParserVCF::parse(r, [&](const Variant &v)
    {
        assert(v.ifs.count("GT"));
        const auto gt = v.ifs.at("GT");
        assert(gt == "HOM" || gt == "HET" || gt == "SOM" || gt == "MSI");
        
        if (gt == "MSI")
        {
            return;
        }
        else if (gt == "HOM")
        {
            x.add(v.name, Mix_1, 1.0);
        }
        else if (gt == "HET")
        {
            x.add(v.name, Mix_1, 0.5);
        }
        else
        {
            assert(!std::isnan(v.allF));
            x.add(v.name, Mix_1, v.allF);
        }
    });

    return x;
}

static Ladder readLadder(const Reader &r, Mixture m, MixtureFormat format, Ladder x = Ladder())
{
    auto parse = [&](const std::string &delim)
    {
        ParserCSV::parse(Reader(r), [&](const ParserCSV::Data &d, Progress i)
        {
            // Don't bother if this is the first line or an invalid line
            if (i == 0 || d.size() <= 1)
            {
                return;
            }

            std::string val;
            
            switch (format)
            {
                case X_M:     { val = d[1]; break; }
                case X_M_X_M: { val = d[3]; break; }
                case M_X_M:   { val = d[2]; break; }
                case X_X_X_M: { val = d[3]; break; }
            }
            
            if (val != "NA")
            {
                x.add(d[0], m, stold(val));
            }
        }, delim);
        
        return x.count();
    };
    
    if (!parse("\t"))
    {
        throw std::runtime_error("No sequin found in the ladder file. Please check and try again.");
    }

    x.src = r.src();
    return x;
}

static Ladder readMix(const Reader &r)
{
    if (countColumns(r) == 2)
    {
        return readLadder(Reader(r), Mix_1, X_M);
    }
    else
    {
        return readLadder(Reader(r), Mix_2, M_X_M, readLadder(Reader(r), Mix_1, X_M));
    }
}

Ladder Standard::readMMix(const Reader &r)
{
    return readMix(r);
}

Ladder Standard::addAF(const Reader &r)
{
    return readAF(r, gen);
}

Ladder Standard::readRGMix(const Reader &r)
{
    auto l1 = readRMix(r);
    Ladder l2;
    
    auto aggregate = [&](Mixture m)
    {
        for (const auto &i : l1.seqs)
        {
            const auto g = RSeq2Std(i);
            
            if (l2.seqs.count(g))
            {
                l2.add(g, m, l1.input(i, m) + l2.input(g, m));
            }
            else
            {
                l2.add(g, m, l1.input(i, m));
            }
        }
    };
    
    aggregate(Mix_1);
    l2.seqs.clear();
    aggregate(Mix_2);
    
    return l2;
}

Ladder Standard::readRLen(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected four columns.");
    return readLadder(Reader(r), Mix_1, X_M);
}

Ladder Standard::readRMix(const Reader &r)
{
    A_CHECK(countColumns(r) == 4, "Invalid mixture file. Expected four columns.");
    auto l = readLadder(Reader(r), Mix_1, M_X_M);
    return readLadder(Reader(r), Mix_2, X_X_X_M, l);
}

template <typename F, typename T = VCFData> T parseVCF2(const Reader &r, F f)
{
    T t;
    
    ParserVCF::parse(r, [&](Variant &x)
    {
        if (f(x))
        {
            t[x.cID].b2v[x.l.start] = x;
            t[x.cID].m2v[x.type()].insert(x);
        }
    });
    
    return t;
}

template <typename F> VCFLadder addVCF(const Reader &r, std::shared_ptr<BedData> rb, F f)
{
    VCFLadder v;
    
    // Only variants fall into the region
    const auto inters = rb ? rb->ginters() : std::map<ChrID, GIntervals<>>();

    v.data = parseVCF2(r, [&](Variant &x)
    {
        A_ASSERT(x.key());
     
        if (rb && (!inters.count(x.cID) || !inters.at(x.cID).overlap(x.l)))
        {
            return false;
        }
        else if (!f(x))
        {
            return false;
        }

        const auto m = std::map<std::string, Genotype>
        {
            { "MSI", Genotype::MSI },
            { "SOM", Genotype::Somatic     },
            { "HOM", Genotype::Homozygous  },
            { "HET", Genotype::Heterzygous },
        };
        
        if (!x.ifs.count("GT") || !m.count(x.ifs.at("GT")))
        {
            throw std::runtime_error(r.src() + " doesn't seem to be a valid VCF reference file. The GT field is not found or invalid.");
        }
        
        assert(!x.name.empty());

        Concent af = 0.0;
        
        switch (x.gt)
        {
            case Genotype::MSI:
            case Genotype::Somatic:     { af = x.allF; break; }
            case Genotype::Homozygous:  { af = 1.0;    break; }
            case Genotype::Heterzygous: { af = 0.5;    break; }
        }
        
        // Update allele frequency ladder
        v.lad.add(x.name, Mix_1, af);
        
        return true;
    });

    return v;
}

VCFLadder Standard::addVCF(const Reader &r, std::shared_ptr<BedData> rb)
{
    auto x = ::addVCF(r, rb, [&](const Variant &)
    {
        return true;
    });
    
    if (x.data.empty())
    {
        throw std::runtime_error("No variant found: " + r.src());
    }
    
    x.src = r.src();
    return x;
}

VCFLadder Standard::addGVCF(const Reader &r, std::shared_ptr<BedData> rb)
{
    auto x = ::addVCF(r, rb, [&](const Variant &v)
    {
        return v.gt != Genotype::Somatic;
    });

    x.src = r.src();
    return x;
}

VCFLadder Standard::addSVCF(const Reader &r, std::shared_ptr<BedData> rb)
{
    auto x =  ::addVCF(r, rb, [&](const Variant &v)
    {
        return v.gt == Genotype::Somatic;
    });

    x.src = r.src();
    return x;
}
