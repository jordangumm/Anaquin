#include "tools/tools.hpp"
#include "data/standard.hpp"
#include "data/resources.hpp"
#include "Genomics/Genomics.hpp"
#include "writers/vcf_writer.hpp"

using namespace Anaquin;

GResource::GResource(const Path &path, const FileName &file, const FileName &ext, HumanAssembly x)
{
    const auto s1 = x == HumanAssembly::chrQS ? "chrQ" : toString(x);
    const auto s2 = toString(x);

    Resource::path = (x == HumanAssembly::None) ? Bundle::latest(path + "/" + file + "_", ext) :
                                                  Bundle::latest(path + "/" + s1 + "/" + file + "_" + s2, ext);
}

LResource::LResource(const Path &path, const FileName &file, const FileName &ext, HumanAssembly x)
{
    Resource::path = Bundle::latest(path + "/" + toString(x) + "/" + file + "_" + toString(x), ext);
}

Bin Anaquin::GBin(const SequinID &x)
{
    auto f = [&](const std::string &s, Bin &r)
    {
             if (s == "IF") { r = Bin::IF; return true; }
        else if (s == "HP") { r = Bin::HP; return true; }
        else if (s == "MS") { r = Bin::MS; return true; }
        else if (s == "MT") { r = Bin::MT; return true; }
        else if (s == "HL") { r = Bin::HL; return true; }
        else if (s == "LD") { r = Bin::LD; return true; }
        else if (s == "SV") { r = Bin::SV; return true; }
        else if (s == "PV") { r = Bin::SV; return true; }
        else if (s == "TL") { r = Bin::SV; return true; }
        else if (s == "IM") { r = Bin::IM; return true; }
        else if (s == "VC") { r = Bin::VC; return true; }
        
        const auto v = Standard::instance().gen.v1()->data.find(x);
        
        if (v)
        {
            switch (v->gt)
            {
                case Genotype::MSI:         { r = Bin::MI; return true; }
                case Genotype::Somatic:     { r = Bin::SO; return true; }
                case Genotype::Homozygous:
                case Genotype::Heterzygous: { r = Bin::GR; return true; }
            }
        }
        
             if (s == "CM") { r = Bin::SO; return true; }
        else if (s == "CV") { r = Bin::SO; return true; }
        else if (s == "CL") { r = Bin::SO; return true; }
        else if (s == "DV") { r = Bin::GR; return true; }
        else if (s == "GV") { r = Bin::GR; return true; }
        
        return false;
    };
    
    Bin r;

    // Consider standard (S0634_CL_047_R) and sequin (S0634_CL_047_R)
    if (f(first(x, "_"), r) || f(second(x, "_"), r))
    {
        return r;
    }

    throw std::runtime_error("Unknown binning for: " + x);
}

bool Anaquin::isHP    (const Bin x) { return x == Bin::HP; }
bool Anaquin::isMS    (const Bin x) { return x == Bin::MS; }
bool Anaquin::isSoma  (const Bin x) { return x == Bin::SO; }
bool Anaquin::isVector(const Bin x) { return x == Bin::VC; }

bool Anaquin::isHP    (const SequinID &x) { return isHP(GBin(x));     }
bool Anaquin::isMS    (const SequinID &x) { return isMS(GBin(x));     }
bool Anaquin::isSoma  (const SequinID &x) { return isSoma(GBin(x));   }
bool Anaquin::isVector(const SequinID &x) { return isVector(GBin(x)); }
bool Anaquin::isGerm  (const SequinID &x) { return isGerm(GBin(x));   }
