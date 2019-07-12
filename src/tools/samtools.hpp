#ifndef SAMTOOLS_HPP
#define SAMTOOLS_HPP

#include <sstream>
#include <htslib/sam.h>

namespace Anaquin
{
    inline std::string bam2qual(const bam1_t *x)
    {
        std::stringstream buf;
        
        for (auto i = 0; i < x->core.l_qseq; ++i)
        {
            buf << (char) (bam_get_qual(x)[i] + 33);
        }
        
        return buf.str();
    }

    inline void bam2qual(char *dst, const bam1_t *x)
    {
        for (auto i = 0; i < x->core.l_qseq; i++, dst++)
        {
            *(dst) = (char) (bam_get_qual(x)[i] + 33);
        }
        
        *(dst) = '\0';
    }

    inline std::string bam2seq(const bam1_t *x)
    {
        std::stringstream buf;

        for (auto i = 0; i < x->core.l_qseq; ++i)
        {
            buf << seq_nt16_str[bam_seqi(bam_get_seq(x), i)];
        }

        return buf.str();
    }

    inline void bam2seq(char *dst, const bam1_t *x)
    {
        for (auto i = 0; i < x->core.l_qseq; i++, dst++)
        {
            *(dst) = seq_nt16_str[bam_seqi(bam_get_seq(x), i)];
        }
        
        *(dst) = '\0';
    }
}

#endif
