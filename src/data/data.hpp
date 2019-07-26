#ifndef DATA_HPP
#define DATA_HPP

#include <map>
#include <cmath>
#include <string>
#include <memory>

namespace Anaquin
{
    typedef unsigned Thread;
    typedef unsigned KmerLen;
    
    typedef std::string Kmer;
    typedef std::string Label;
    typedef std::string Sequence;
    
    typedef long long KMers;
    typedef long long Reads;
    
    typedef double Coverage;
    typedef double Measured;
    typedef double KMCoverage;
    
    typedef double Proportion;
    typedef long double Concent;
    typedef long double Probability;
    
    typedef std::string Name;
    typedef std::string ChrID;
    typedef std::string SequinID;
    typedef std::string StandardID;
    
    typedef std::string Path;
    typedef std::string Label;
    typedef std::string Token;
    typedef std::string Units;
    typedef std::string Command;
    typedef std::string Scripts;
    typedef std::string Version;
    typedef std::string FileName;
    typedef std::string ReadName;
    
    typedef long long Base;
    typedef long long Depth;
    typedef long long Counts;
    
    const auto NO_CALIBRATION = -1.0;

    enum Variation
    {
        SNP,
        Insertion,
        Deletion,
        Inversion,
        Duplication,
    };
    
    enum Genotype
    {
        Somatic,
        MSI,
        Homozygous,
        Heterzygous
    };
    
    enum class Mutation
    {
        Germline,
        Somatic
    };

    enum HumanAssembly
    {
        gr37  = 1,
        gr38  = 2,
        hg19  = 3,
        hg38  = 4,
        chrQS = 5,
        hseq  = 6,
        None  = 7
    };
    
    class Confusion
    {
        public:
        
            void operator+=(const Confusion &m)
            {
                _tp += m.tp();
                _fp += m.fp();
                _fn += m.fn();
                _nr += m.nr();
                _nq += m.nq();
            }
        
            inline Counts &tp() const { return _tp; }
            inline Counts &fp() const { return _fp; }
            inline Counts &fn() const { return _fn; }
            inline Counts &nr() const { return _nr; }
            inline Counts &nq() const { return _nq; }
        
            // Sensitivity (true positive rate)
            inline Proportion sn() const
            {
                return (_tp + _fn) ? static_cast<Proportion>(_tp) / (_tp + _fn) : NAN;
            }
        
            // Precision (positive predictive value)
            inline Proportion pc() const
            {
                return (_tp + _fp) ? static_cast<Proportion>(_tp) / (_tp + _fp) : NAN;
            }
        
            // False discovery (false discovery rate)
            inline Proportion fdr() const
            {
                return 1.0 - pc();
            }
        
            inline Proportion F1() const
            {
                return 2.0 * ((pc() * sn()) / (pc() + sn()));
            }
        
        private:
        
            mutable Counts _fp = 0;
            mutable Counts _tp = 0;
            mutable Counts _fn = 0;
            mutable Counts _nq = 0;
            mutable Counts _nr = 0;
    };
    
    enum Bin
    {
        IF, // Info
        MT, // Mitro
        MS, // Micro
        HL,
        HP,
        LD,
        SV,
        IM, // Immune
        ES,
        GR, // Germline (or just sequins)
        SO, // Somatic
        MI, // MSI
        VC
    };

    template <typename F> void forBin(F f)
    {
        for (auto i = 0; i <= (int) VC; i++)
        {
            f((Bin) i);
        }
    }

    enum class Product
    {
        Genomics,
        RNA,
        Meta,
    };
    
    struct AlignedRead
    {
        void *b; // bam1_t
        bool rc; // Reverse complement?
    };
    
    typedef std::string Line;
    typedef long long Progress;
    typedef std::pair<std::shared_ptr<AlignedRead>, std::shared_ptr<AlignedRead>> PairedRead;
    typedef std::map<ReadName, PairedRead> AlignedReads;

    /*
     * Represents a mathced element that can be identified
     */
    
    struct Matched
    {
        virtual std::string name() const = 0;
    };
    
    enum Mixture
    {
        Mix_1,
        Mix_2,
        Mix_3
    };

    struct Limit
    {
        SequinID id;
        
        // Expected concentration for the sequin
        Concent abund = NAN;
    };
    
    struct KMInfo
    {
        Kmer kmer;
        Counts abund = 0;
    };
    
    const std::string MISSING = "NA";
}

#endif
