#ifndef RNA_HPP
#define RNA_HPP

#include "data/bundle.hpp"

namespace Anaquin
{
    #define DEFAUL_RNA_CALIB_FRAC 0.01

    struct RResource : public Resource
    {
        RResource(const Path &path, const FileName &file, const FileName &ext)
        {
            Resource::path = Bundle::latest(path + "/" + file, ext);
        }
    };

    inline Resource RNAMix(const Path &p)
    {
        return RResource(p + "/transcriptome", "rnasequin_isoforms_", ".tsv");
    }
    
    inline Resource RNAFA(const Path &p)
    {
        return RResource(p + "/transcriptome", "rnasequin_sequences_", ".fa");
    }

    inline Bin RBin(const SequinID &) { return GR; }
    
    inline bool RValid(Bin x) { return x == GR || x == ES; }
    
    inline StandardID RSeq2Std(const SequinID &x) { return noLast(x, "_"); }
}

#endif
