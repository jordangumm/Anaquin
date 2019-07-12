#ifndef META_HPP
#define META_HPP

#include "data/bundle.hpp"

namespace Anaquin
{
    #define DEFAUL_META_CALIB_FRAC 0.01

    struct MResource : public Resource
    {
        MResource(const Path &path, const FileName &file, const FileName &ext)
        {
            Resource::path = Bundle::latest(path + "/" + file, ext);
        }
    };
    
    inline Resource MetaFA(const Path &p)
    {
        return MResource(p + "/metagenome", "metasequin_sequences_", ".fa");
    }

    inline Resource MetaMix(const Path &p)
    {
        return MResource(p + "/metagenome", "metasequin_abundance_", ".tsv");
    }
    
    inline Bin MBin(const SequinID &x)
    {
        if (second(x, "_") == "VC")
        {
            return VC;
        }
        else if (first(x, "_") == "LD" || second(x, "_") == "LD")
        {
            return LD;
        }
        else if ((first(x, "_") == "IF" || second(x, "_") == "IF") ||
                 (first(x, "_") == "UQ" || second(x, "_") == "UQ"))
        {
            return IF;
        }
        
        return GR;
    }
        
    inline bool MValid(Bin x) { return x == GR || x == ES || x == LD || x == IF || x == VC; }

    inline StandardID MSeq2Std(const SequinID &x) { return isSubstr(x, "LD_") ? noLast(noFirst(x, "_"), "_") : x; }
}

#endif
