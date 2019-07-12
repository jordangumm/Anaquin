#ifndef GENOMICS_HPP
#define GENOMICS_HPP

#include "data/data.hpp"
#include "data/bundle.hpp"
#include "tools/tools.hpp"

namespace Anaquin
{
    struct GResource : public Resource
    {
        GResource(const Path &, const FileName &, const FileName &, HumanAssembly);
    };

    struct LResource : public Resource
    {
        LResource(const Path &, const FileName &, const FileName &, HumanAssembly);
    };

    inline Resource GInfoCode(const Path &p)
    {
        return LResource(p + "/synthetic", "sequin_barcodes", ".tsv", HumanAssembly::None);
    }

    inline Resource GSeqFA(const Path &p)
    {
        return GResource(p + "/genome", "sequin_sequences", ".fa", HumanAssembly::None);
    }

    inline Resource GRegionBED(const Path &p, HumanAssembly x)
    {
        return GResource(p + "/genome", "sequin_regions", ".bed", x);
    }

    inline Resource GVarVCF(const Path &p, HumanAssembly x)
    {
        return GResource(p + "/genome", "sequin_smallvariants", ".vcf", x);
    }

    inline Resource GFeatBED(const Path &p, HumanAssembly x)
    {
        return GResource(p + "/genome", "sequin_features", ".bed", x);
    }

    inline Resource GAttrBED(const Path &p)
    {
        return GResource(p + "/genome", "sequin_attributes", ".tsv", HumanAssembly::None);
    }

    inline Resource GSynTSV(const Path &p)
    {
        return GResource(p + "/synthetic", "synthetic_ladder", ".tsv", HumanAssembly::None);
    }

    Bin GBin(const SequinID &);

    inline std::string gt2str(Genotype x)
    {
        switch (x)
        {
            case Genotype::MSI:         { return "MSI";          }
            case Genotype::Somatic:     { return "Somatic";      }
            case Genotype::Homozygous:  { return "Homozygous";   }
            case Genotype::Heterzygous: { return "Heterozygous"; }
        }
    }

    inline std::string var2str(Variation x)
    {
        switch (x)
        {
            case Variation::SNP:         { return "SNP";         }
            case Variation::Deletion:    { return "Deletion";    }
            case Variation::Insertion:   { return "Insertion";   }
            case Variation::Duplication: { return "Duplication"; }
            case Variation::Inversion:   { return "Inversion";   }
        }
    }

    inline StandardID GSeq2Std(const SequinID &x)
    {
        return noLast(noFirst(x, "_"), "_");
    }
    
    inline bool GValid(Bin) { return true; }

    inline SequinID __hack__(const SequinID &x) { return remove(x, "_A") + "_01"; }
    
    bool isHP(const Bin x);
    bool isHP(const SequinID &x);
    
    bool isMS(const Bin x);
    bool isMS(const SequinID &x);

    bool isSoma(const Bin x);
    bool isSoma(const SequinID &x);

    bool isVector(const Bin x);
    bool isVector(const SequinID &x);

    bool isGerm(const SequinID &x);
    inline bool isGerm(Bin x) { return x == Bin::GR;  }
}

#endif
