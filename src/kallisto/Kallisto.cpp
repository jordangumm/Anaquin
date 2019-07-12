#include <assert.h>
#include <unistd.h>
#include <iostream>
#include <sys/types.h>
#include <htslib/bgzf.h>
#include "RNA/RNA.hpp"
#include "KmerIndex.h"
#include "Kallisto.hpp"
#include "Meta/Meta.hpp"
#include "data/fastq.hpp"
#include "ProcessReads.h"
#include "tools/tools.hpp"
#include "KmerIterator.hpp"
#include "tools/errors.hpp"
#include <boost/format.hpp>
#include "tools/samtools.hpp"
#include "Genomics/Genomics.hpp"
#include "writers/bam_writer.hpp"

// Combined?
extern int __Combined__;

// Defined in ProcessReads.cpp
extern bam_hdr_t *__CombinedBAMHeader__;

using namespace Anaquin;

// Cache for combined BAM
static std::map<Thread, std::unique_ptr<char>> s1_, s2_, q1_, q2_;

// Extremly unlikely but it's only 1MB
static const int MAX_READ_LENGTH = 1000000;

class Partition
{
    public:
        Partition(const Path &, Thread);
        ~Partition() { close(); }

        void close()
        {
            for (auto &i : o1)   { i.second->close();    }
            for (auto &i : o2)   { i.second->close();    }
            for (auto &i : bamW) { bgzf_close(i.second); }
            bamW.clear();
        }
    
        void writeB(Bin k, const bam1_t *x)
        {
            bam_write1(bamW.at(k), x);
        }
    
        void write(Bin, const char *, const char *, const char *, void *,
                        const char *, const char *, const char *, void *, bool);
    
        /*
         * Partitioned reads to FASTQ
         */
    
        std::map<Bin, FileName> c1, c2;
        std::map<Bin, std::shared_ptr<std::ofstream>> o1, o2;

        /*
         * Partitioned reads to BAM
         */
    
        std::map<Bin, BGZF *>   bamW;
        std::map<Bin, FileName> bamF;
};

struct Index
{
    // Partition for this index
    std::map<Thread, std::shared_ptr<Partition>> part;
    
    // Running statistics by threads
    std::map<Thread, KStats> stats;
    
    // Kallisto index
    std::shared_ptr<KmerIndex> index;
};

typedef std::map<unsigned, std::shared_ptr<Partition>> PartitionReads;

static Index __R__;

// Temp directory for the current exectution
static Path __tmp__;

static KOptions __o__;

// Describe the FASTQ file(s)
FASTQ __KFQInfo__;

void Partition::write(Bin k, const char *r1, const char *s1, const char *q1, void *b1,
                             const char *r2, const char *s2, const char *q2, void *b2, bool rev)
{
    auto writeF = [&](const char *r, const char *s, const char *q, std::shared_ptr<std::ofstream> o, bool rev)
    {
        assert(o);
        std::string ss = s;
        if (rev) { complement(ss); }
        
        std::stringstream str;
        str << "@" + std::string(r) + "\n";
        str << ss + "\n";
        str << "+\n";
        str << std::string(q) + "\n";
        
        const auto x = compressGZ(str.str());
        o->write(x.data(), x.size());
    };

    if (__o__.writeBAM())
    {
        if (!__o__.onlySeqLad || k == GR || k == LD)
        {
            writeB(k, (bam1_t *) b1);
            writeB(k, (bam1_t *) b2);
        }
    }
    else
    {
        writeF(r1, s1, q1, o1[k], rev);
        writeF(r2, s2, q2, o2[k], rev);
    }
}

Partition::Partition(const Path &p, Thread tID)
{
    forBin([&](Bin c)
    {
        bool x;
        
        switch (__o__.prod)
        {
            case Product::RNA:      { x = RValid(c); break; }
            case Product::Meta:     { x = MValid(c); break; }
            case Product::Genomics: { x = GValid(c); break; }
        }
        
        if (__o__.writeBAM())
        {
            assert(__CombinedBAMHeader__);
            
            // Never write sample reads for BAM (because SCombine() did it)
            if (c != ES)
            {
                if (!__o__.onlySeqLad || c == GR || c == LD)
                {
                    bamF[c] = p + "/partition_B" + std::to_string(c) + "_" + toString(tID) + ".bam.tmp";
                    bam_hdr_write(bamW[c] = bgzf_open(bamF[c].c_str(), "w"), __CombinedBAMHeader__);
                }
            }
        }
        else if (x)
        {
            const auto f1 = p + "/partition_B" + std::to_string(c) + "_" + toString(tID) + "_1.fq.gz.tmp";
            const auto f2 = p + "/partition_B" + std::to_string(c) + "_" + toString(tID) + "_2.fq.gz.tmp";
            o1[c] = std::shared_ptr<std::ofstream>(new std::ofstream(c1[c] = f1, std::ios::binary | std::ios::out));
            o2[c] = std::shared_ptr<std::ofstream>(new std::ofstream(c2[c] = f2, std::ios::binary | std::ios::out));
        }
    });
}

static std::shared_ptr<KmerIndex> KIndex(const FileName &file, KmerLen k)
{
    ProgramOptions o;
    o.k = k;
    o.index = file;
    
    auto x = std::shared_ptr<KmerIndex>(new KmerIndex(o));
    x->load(o);
    return x;
}

// Build Kallisto index from a FASTA file
static FileName KBuildIndex(const FileName &file, unsigned k)
{
    ProgramOptions opt;
    
    opt.k = ::Kmer::k = k;
    opt.index = tmpFile();
    opt.transfasta.push_back(file);
    
    KmerIndex index(opt);
    index.BuildTranscripts(opt);
    index.write(opt.index);
    
    if (index.dbGraph.contigs.empty() || !index.kmap.size())
    {
        throw std::runtime_error("Failed to build index for " + file);
    }
    
    return opt.index;
}

// Initialise Kallisto. Required before running anything else
static void KInit(const KOptions &o)
{
    __KFQInfo__ = FASTQ();
    
    // Always initialize in case it has been used (like ladder recalibration)
    __R__ = Index();

    __o__ = o;
    assert(__o__.prod == Product::RNA || __o__.prod == Product::Meta || __o__.prod == Product::Genomics);

    // Running combined?
    __Combined__ = o.bam ? 1 : 0;

    __tmp__ = tmpPath();
    removeD(__tmp__);
    createD(__tmp__);

    __R__.index = KIndex(KBuildIndex(o.index, o.k), o.k);

    if (__R__.index->target_names_.empty())
    {
        throw std::runtime_error("No sequence found in reference index: " + o.index);
    }

    if (o.bam)
    {
        for (auto i = 0; i < o.thr; i++)
        {
            s1_[i] = std::unique_ptr<char>(new char[MAX_READ_LENGTH]);
            s2_[i] = std::unique_ptr<char>(new char[MAX_READ_LENGTH]);
            q1_[i] = std::unique_ptr<char>(new char[MAX_READ_LENGTH]);
            q2_[i] = std::unique_ptr<char>(new char[MAX_READ_LENGTH]);
        }
    }

    for (auto i = 0u; i < __R__.index->target_names_.size(); i++)
    {
        const auto &seq = __R__.index->target_names_[i];
        
        StandardID std;
        
        switch (__o__.prod)
        {
            case Product::RNA:      { std = RSeq2Std(seq); break; }
            case Product::Meta:     { std = MSeq2Std(seq); break; }
            case Product::Genomics: { std = GSeq2Std(seq); break; }
        }

        assert(!std.empty());

        for (auto j = 0; j < o.thr; j++)
        {
            __R__.stats[j].stds.insert(std);
            __R__.stats[j].seqs.insert(seq);
            
            if (__o__.prod == Product::Genomics)
            {
                auto isRSeq = [&](const SequinID &x)
                {
                    return isEnd(x, "_R");
                };
                
                if (isRSeq(seq))
                {
                    __R__.stats[j].rSeqs[std] = seq;
                }
                else
                {
                    __R__.stats[j].aSeqs[std].insert(seq);
                }
            }
        }
    }
}

struct EdgeMatch
{
    // n-th k-mer
    Base i;
    
    // Position matching in the reference
    Base p;
};

static bool matchIndex(Counts nMatch, Counts nNMatch)
{
    const auto mp = static_cast<float>(nMatch) / (nMatch + nNMatch);
    return nMatch > nNMatch || mp >= __o__.rule;
}

// Match a read by indexing
static Bin KMatch(Bin (*f)(const std::string &),
                  Thread tID,
                  const char *,
                  const char *s,
                  void *b,
                  std::map<SequinID, Counts> &ms,
                  bool &isForw,
                  Base &start,
                  Base &end)
{
    // Our only index
    auto &x = __R__;

    if (!x.part.count(tID))
    {
        __R__.part.insert(std::pair<unsigned, std::shared_ptr<Partition>>(
                          tID, std::shared_ptr<Partition>(
                          new Partition(__tmp__, tID))));
    }
    
    std::string tmp;
    
    if (__o__.flipBefore)
    {
        tmp = s;
        complement(tmp);
        s = tmp.data();
    }

    // Number of k-mers in the read
    const auto n = strlen(s) - __o__.k + 1;

    auto __match__ = [&](std::shared_ptr<KmerIndex> index)
    {
        ms.clear();
        KmerIterator ki(s), ke;
        
        Counts nMatch  = 0;
        Counts nNMatch = 0;
        
        SequinKM shared, uniqs;
        
        auto j = 0;
        bool isLastMatched = false;
        
        // Edge matching for first and last
        std::shared_ptr<EdgeMatch> fm, lm;
        
        for (int i = 0; ki != ke; ++i, ++ki)
        {
            // Never skip the first k-mer
            if (i && j != __o__.skipKM)
            {
                j++;
                continue;
            }
            
            j = 0;
            auto search = index->kmap.find(ki->first.rep());

            // Where the k-mer starts
            const auto l = i;
            
            // Where the k-mer finishes
            //const auto u = l + __o__.k - 1;

            // Can we find this k-mer?
            if (search != index->kmap.end())
            {
                isLastMatched = true;
                
                nMatch++;
                const auto km = search->first.toString();

                std::set<SequinID> seqs;
                for (const auto &tran : index->dbGraph.contigs[search->second.contig].transcripts)
                {
                    seqs.insert(index->target_names_[tran.trid]);
                    
                    auto findPos = [&]()
                    {
                        auto x = index->findPosition(tran.trid, search->first, search->second);
                        return x.second ? x.first : x.first - __o__.k;
                    };
                    
                    if (!fm)
                    {
                        fm = std::shared_ptr<EdgeMatch>(new EdgeMatch());
                        fm->i = i; fm->p = findPos();
                    }
                    else
                    {
                        if (!lm) { lm = std::shared_ptr<EdgeMatch>(new EdgeMatch()); }
                        lm->i = i; lm->p = findPos();
                    }
                }

                // Unique matching?
                const auto isUniq = seqs.size() == 1;
                
                for (const auto &seq : seqs)
                {
                    ms[seq]++;

                    if (!isUniq) { shared[seq][km]++; }
                    else         { uniqs [seq][km]++;  }
                }
            }
            else
            {
                isLastMatched = false;
                nNMatch++;
            }
        }

        auto copy = [&](const SequinKM &x1, SequinKM &x2)
        {
            for (const auto &i : x1)
            {
                for (const auto &j : i.second)
                {
                    x2[i.first][j.first] += j.second;
                }
            }
        };
        
        // Only matched a single k-mer?
        if (fm && !lm) { lm = fm; }

        // Match this index?
        if (matchIndex(nMatch, nNMatch))
        {
            copy(uniqs,  x.stats[tID].uniqs);
            copy(shared, x.stats[tID].shared);
            
            x.stats[tID].nMatch++;
            
            // How many matching k-mers?
            x.stats[tID].nMKMatch += nMatch;
            
            // How many non-matching k-mers?
            x.stats[tID].nNMKMatch += nNMatch;

            if (!fm || !lm)
            {
                // Assume the whole sequence
                start = 0;
                
                // Assume the whole sequence
                end = strlen(s);
            }
            else
            {
                lm = !lm ? fm : lm;
                
                // Forward strand?
                isForw = fm->p < lm->p;
                
                // Estimated beginning
                start = isForw ? (fm->p - fm->i) : lm->p - (n - lm->i - 1);
                
                // Estimated ending
                end = isForw ? (lm->p + __o__.k + (n - lm->i - 1)) : fm->p + (__o__.k - 1) + (n - lm->i - 1);
            }

            return true;
        }
        
        // Doesn't match this index
        else
        {
            ms.clear();
            x.stats[tID].nNMatch++;
            return false;
        }
    };
    
    // Force reads if we're running combined mode
    return (__match__(x.index)) ? (ms.empty() ? GR: f(max(ms))) : (__Combined__ == 1 ? GR : ES);
}

void KPartition(Thread tID,
                const char *r1,
                const char *s1,
                const char *q1,
                void       *b1,
                bool       rc1,
                const char *r2,
                const char *s2,
                const char *q2,
                void       *b2,
                bool       rc2)
{
    if ((tID == 0) && !__KFQInfo__.heads())
    {
        __KFQInfo__.addHead(r1);
    }
    
    std::string tmp1, tmp2, tmp3, tmp4;
    
    if (__o__.bam)
    {
        assert(!s1 && !q1 && !s2 && !q2);
        assert(s1_.count(tID) && q1_.count(tID));
        assert(s2_.count(tID) && q2_.count(tID));

        bam2seq (s1_.at(tID).get(), static_cast<const bam1_t *>(b1));
        bam2qual(q1_.at(tID).get(), static_cast<const bam1_t *>(b1));
        bam2seq (s2_.at(tID).get(), static_cast<const bam1_t *>(b2));
        bam2qual(q2_.at(tID).get(), static_cast<const bam1_t *>(b2));
        
        s1 = s1_[tID].get(); q1 = q1_[tID].get();
        s2 = s2_[tID].get(); q2 = q2_[tID].get();

        if (rc1)
        {
            tmp1 = std::string(q1);
            std::reverse(tmp1.begin(), tmp1.end());
            q1 = tmp1.c_str();
            tmp2 = s1;
            tmp2 = revcomp(tmp2);
            s1 = tmp2.c_str();
        }
        
        if (rc2)
        {
            tmp3 = std::string(q2);
            std::reverse(tmp3.begin(), tmp3.end());
            q2 = tmp3.c_str();
            tmp4 = s2;
            tmp4 = revcomp(tmp4);
            s2 = tmp4.c_str();
        }

        assert(strlen(s1) && strlen(s2));
        assert(strlen(s1) == strlen(q1));
        assert(strlen(s2) == strlen(q2));
    }

    // Counting for all matches (all k-mers)
    std::map<SequinID, Counts> R1, R2, F1, F2;

    // Forward orientation?
    bool isF1, isF2;
    
    // Position for trimming
    Base start1, end1, start2, end2;
    
    auto match = [&](Index &i, std::map<SequinID, Counts> &m1, std::map<SequinID, Counts> &m2)
    {
        switch (__o__.prod)
        {
            case Product::Genomics:
            {
                return std::pair<Bin, Bin>(KMatch(GBin, tID, r1, s1, b1, m1, isF1, start1, end1),
                                           KMatch(GBin, tID, r2, s2, b2, m2, isF2, start2, end2));
            }
                
            case Product::RNA:
            {
                return std::pair<Bin, Bin>(KMatch(RBin, tID, r1, s1, b1, m1, isF1, start1, end1),
                                           KMatch(RBin, tID, r2, s2, b2, m2, isF2, start2, end2));
            }
                
            case Product::Meta:
            {
                return std::pair<Bin, Bin>(KMatch(MBin, tID, r1, s1, b1, m1, isF1, start1, end1),
                                           KMatch(MBin, tID, r2, s2, b2, m2, isF2, start2, end2));
            }
        }
    };

    auto fill = [&](Bin k, Index &i, const std::map<SequinID, Counts> &m1, const std::map<SequinID, Counts> &m2)
    {
        const auto x1 = !m1.empty() ? max(m1) : "-";
        const auto x2 = !m2.empty() ? max(m2) : "-";
        
        auto addRead = [&]()
        {
            auto r1_ = first(r1, " ");
                 r1_ = isSubstr(r1_, "/") ? noLast(r1_, "/") : r1_;
            auto r2_ = first(r2, " ");
                 r2_ = isSubstr(r2_, "/") ? noLast(r2_, "/") : r2_;

            assert(!r1_.empty() && !r2_.empty());
            
            i.stats[tID].r1[k][r1_] = x1;
            i.stats[tID].r2[k][r2_] = x2;
        };
        
        // Always required for ladder calibration
        if (k == Bin::LD)
        {
            addRead();
        }

        i.stats[tID].sqc[x1]++;
        i.stats[tID].sqc[x2]++;
        i.stats[tID].c1[k]++;
        i.stats[tID].c2[k]++;

        auto flip = !(k == ES || k == VC || k == IF || k == LD);
       
        if (__o__.prod == Product::RNA)
        {
           flip = false;
        }
 
        if (!__o__.flip)
        {
            flip = false;
        }

        if (i.part.count(tID))
        {
            i.part[tID]->write(k, r1, s1, q1, b1, r2, s2, q2, b2, flip);
        }
    };

    auto write = [&](Bin x1, Bin x2)
    {
        if (__Combined__ == 0 && ((x1 != Bin::SV && x2 == Bin::SV) || (x1 == Bin::SV && x2 != Bin::SV)))
        {
            fill(ES, __R__, R1, R2);
        }
        else if (x1 == Bin::VC || x2 == Bin::VC)
        {
            fill(VC, __R__, R1, R2);
        }
        else
        {
            if (__Combined__ == 0 && x1 != x2)
            {
                // Write ambigious reads to the bin for sample reads
                fill(ES, __R__, R1, R2);
            }
            else
            {
                fill(x1, __R__, R1, R2);
            }
        }
    };
    
    const auto R = match(__R__, R1, R2);
    write(R.first, R.second);
}

// Merge the internal Kallisto's results to more accessible KStats
static KStats KMerge()
{
    auto &x = __R__;
    auto r = x.stats[0];
    
    // Close all connections
    for (const auto &i : __R__.part) { i.second->close(); }

    // For each thread, merge with all previous threads...
    for (auto i = 1u; i < x.stats.size(); i++)
    {
        auto addX1X2 = [](KStats &x1, KStats &x2)
        {
            assert(x1.aSeqs.size() == x2.aSeqs.size());
            
            // We'll just need to add x2
            KStats x = x1;
            
            auto f = [&](KStats &x1, KStats &x2)
            {
                x1.nMatch    += x2.nMatch;
                x1.nNMatch   += x2.nNMatch;
                x1.nMKMatch  += x2.nMKMatch;
                x1.nNMKMatch += x2.nNMKMatch;
                
                x1.uniqs  = add(x1.uniqs, x2.uniqs);
                x1.shared = add(x1.shared, x2.shared);
                
                x1.c1  = add(x1.c1,  x2.c1);
                x1.c2  = add(x1.c2,  x2.c2);
                x1.sqc = add(x1.sqc, x2.sqc);
                
                for (auto i = 0; i <= (int)VC; i++)
                {
                    x1.r1[(Bin)i] = add(x1.r1[(Bin)i], x2.r1[(Bin)i]);
                    x1.r2[(Bin)i] = add(x1.r2[(Bin)i], x2.r2[(Bin)i]);
                }
            };
            
            f(x, x2);
            return x;
        };
        
        r = addX1X2(r, x.stats[i]);
    }
    
    for (const auto &i : x.part)
    {
        forBin([&](Bin c)
        {
            if (__o__.writeBAM())
            {
                if (c != ES)
                {
                    if (!__o__.onlySeqLad || c == GR || c == LD)
                    {
                        r.f1[c].push_back(i.second->bamF.at(c));
                    }
                }
            }
            else
            {
                if (r.c1.count(c) && r.c2.count(c))
                {
                    r.f1[c].push_back(i.second->c1.at(c));
                    r.f2[c].push_back(i.second->c2.at(c));
                }
            }
        });
    }
    
    s1_.clear(); s2_.clear(); q1_.clear(); q2_.clear();
    r.work = __tmp__; __tmp__.clear();
    return r;
};

KStats Kallisto(const FileName &f1, const FileName &f2, const KOptions &o)
{
    KInit(o);
    ProgramOptions opt;
    
    opt.k = ::Kmer::k = o.k;
    opt.index = o.index;
    opt.threads = o.thr;
    opt.files.push_back(f1);
    opt.files.push_back(f2);

    // Required for feteching read names
    opt.fusion = true;
    
    KmerIndex index(opt);
    MinCollector collection(index, opt);
    ProcessReads(index, opt, collection);
    
    return KMerge();
}
