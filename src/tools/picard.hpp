#ifndef PICARD_HPP
#define PICARD_HPP

#include <map>
#include <htslib/sam.h>
#include "data/standard.hpp"
#include "data/alignment.hpp"

namespace Anaquin
{
    enum class SNPBin
    {
        AC, AT, AG,
        CA, CT, CG,
        TA, TG, TC,
        GA, GC, GT,
        Match
    };
    
    struct PicardOption
    {
        bool ignoreSkipClip = false;
    };
    
    struct Picard
    {
        Picard(const std::map<ChrID, Sequence> &seqs) : seqs(seqs)
        {
            for (const auto &i : seqs)
            {
                snps[i.first]; ins[i.first]; dls[i.first];
            }
        }
        
        void analyze(const ChrID &cID, const Alignment &x, const PicardOption &o = PicardOption())
        {
            assert(seqs.count(cID));
            assert(!x.seq.empty());
            
            const auto seq = seqs.at(cID);
            
            const auto &r = Standard::instance().gen;
            assert(r.v4());
            
            // Relative to the query sequence
            Base i = 0;
            
            // Relative to reference genome
            Base j = x.l.start-1;
            
            for (auto c: x.cigars)
            {
                if (o.ignoreSkipClip && (c.first == BAM_CREF_SKIP || c.first == BAM_CSOFT_CLIP || c.first == BAM_CHARD_CLIP))
                {
                    return;
                }
            }
            
            for (auto c: x.cigars)
            {
                const auto val = c.second;
                
                switch (c.first)
                {
                    case BAM_CMATCH:
                    {
                        for (auto k = 0; k < val; k++)
                        {
                            assert(j < seq.size());
                            
                            if (seq.at(j) == x.seq[i])
                            {
                                snps[seq][j][SNPBin::Match]++;
                            }
                            else
                            {
                                assert(i < x.seq.size());
                                
                                if (r.v4()->data.findVar(x.cID, Locus(j+1, j+1)))
                                {
                                    // Skip variant allele
                                    continue;
                                }
                                
                                const auto r = toupper(seq.at(j));
                                const auto q = toupper(x.seq[i]);
                                
                                if (r != 'N' && q != 'N')
                                {
                                    assert(r == 'A' || r == 'C' || r == 'T' || r == 'G');
                                    assert(q == 'A' || q == 'C' || q == 'T' || q == 'G');
                                    
                                    SNPBin c = SNPBin::AC;
                                    
                                    if (r == 'A' && q == 'C') { c = SNPBin::AC; }
                                    if (r == 'A' && q == 'T') { c = SNPBin::AT; }
                                    if (r == 'A' && q == 'G') { c = SNPBin::AG; }
                                    if (r == 'C' && q == 'A') { c = SNPBin::CA; }
                                    if (r == 'C' && q == 'T') { c = SNPBin::CT; }
                                    if (r == 'C' && q == 'G') { c = SNPBin::CG; }
                                    if (r == 'T' && q == 'A') { c = SNPBin::TA; }
                                    if (r == 'T' && q == 'G') { c = SNPBin::TG; }
                                    if (r == 'T' && q == 'C') { c = SNPBin::TC; }
                                    if (r == 'G' && q == 'A') { c = SNPBin::GA; }
                                    if (r == 'G' && q == 'C') { c = SNPBin::GC; }
                                    if (r == 'G' && q == 'T') { c = SNPBin::GT; }
                                    
                                    snps[seq][j][c]++;
                                }
                            }
                            
                            i++;
                            j++;
                        }
                        
                        break;
                    }
                        
                    case BAM_CINS: { ins[seq][j][val]++; i += val; break; }
                    case BAM_CDEL: { dls[seq][j][val]++; j += val; break; }
                        
                    case BAM_CREF_SKIP:
                    case BAM_CSOFT_CLIP: { skips += val; i += val; break; }
                    case BAM_CHARD_CLIP: { skips += val; break; }
                        
                    case BAM_CPAD:   { assert(false); break; }
                    case BAM_CDIFF:  { assert(false); break; }
                    case BAM_CBACK:  { assert(false); break; }
                    case BAM_CEQUAL: { assert(false); break; }
                        
                    default:
                    {
                        throw std::runtime_error("Unknown " + std::to_string(c.first));
                    }
                }
            }
        }
        
        std::string report() const;
        
        typedef Base RBase;
        typedef Base Length;
        
        typedef std::map<ChrID, std::map<RBase, std::map<SNPBin, Counts>>> SNPData;
        typedef std::map<ChrID, std::map<RBase, std::map<Length, Counts>>> InsDelData;
        
        mutable SNPData snps;
        mutable InsDelData ins, dls;
        
        inline Counts totalS()
        {
            Counts n = 0;
            
            for (const auto &i : snps)
            {
                for (const auto &j : i.second)
                {
                    for (const auto &k : j.second) { n += k.second; }
                }
            }
            
            return n;
        }
        
        // Sums all counts for a SNP bin
        inline Counts sumSB(SNPBin b)
        {
            Counts n = 0;
            
            for (const auto &i : snps)
            {
                for (const auto &j : i.second)
                {
                    for (const auto &k : j.second)
                    {
                        if (k.first == b)
                        {
                            n += k.second;
                        }
                    }
                }
            }
            
            return n;
        }

        // Sums all counts for either insertion or deletion
        static Counts sumL(const InsDelData &x, Length l)
        {
            Counts n = 0;
            
            for (const auto &i : x)
            {
                for (const auto &j : i.second)
                {
                    for (const auto &k : j.second)
                    {
                        if (l == 0 || k.first == l)
                        {
                            n += k.second;
                        }
                    }
                }
            }
            
            return n;
        }
    
        // Sum of deletions (0 for all sizes)
        inline Counts sumDL(Length l = 0) { return Picard::sumL(dls, l); }
        
        // Sum of insertions (0 for all sizes)
        inline Counts sumIL(Length l = 0) { return Picard::sumL(ins, l); }

        Counts skips = 0;

        // Reference sequences
        const std::map<ChrID, Sequence> seqs;
    };
}

#endif
