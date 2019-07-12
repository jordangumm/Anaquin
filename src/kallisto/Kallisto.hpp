#ifndef KALLISTO_HPP
#define KALLISTO_HPP

#include <map>
#include <memory>
#include "tools/tools.hpp"

namespace Anaquin
{
    typedef std::map<SequinID, std::map<Kmer, Counts>> SequinKM;

    // Default k-mer length
    const KmerLen K_DEFAULT_L = 23;
    
    // Default classification rule
    const Proportion K_DEFAULT_R  = 0.20;

    /*
     * Statistics from Kallisto
     */
    
    struct KStats
    {
        Path work;
        
        // Counts for the categories
        std::map<Bin, Counts> c1, c2;
        
        // Output file names (f2 might not be defined)
        std::map<Bin, std::vector<FileName>> f1, f2;
 
        // Total number of reads
        inline Counts total() const { return sum(c1); }

        // Number of reads for a bin
        inline Counts binN(Bin x) const { return c1.count(x) ? c1.at(x) : 0; }

        // Percentage of reads for a bin
        inline Proportion binP(Bin x) const { return c1.count(x) ? (c1.at(x) / (float) total()) : 0.0; }
        
        std::set<StandardID> stds, seqs;
        std::map<StandardID, SequinID> rSeqs;
        std::map<StandardID, std::set<SequinID>> aSeqs;
        
        // Measured counts for unique k-mers
        SequinKM uniqs;
        
        // Measured counts for shared k-mers
        SequinKM shared;
        
        // Number of reads matching this index (e.g. sequins)
        Counts nMatch = 0;
        
        // Number of reads not matching this index (e.g. genome)
        Counts nNMatch = 0;
        
        // Number of matching k-mers in all reads matching this index
        Counts nMKMatch = 0;
        
        // Number of non-matching k-mers in all reads matching this index
        Counts nNMKMatch = 0;
        
        // Number of reads for each sequin (2x for paired-end)
        std::map<SequinID, Counts> sqc;
        
        // Information about paired-end reads
        std::map<Bin, std::map<ReadName, SequinID>> r1, r2;
    };

    struct KOptions
    {
        KOptions() : flipBefore(false), flip(true), onlySeqLad(false), forceFQ(false) {}
        
        inline bool writeBAM() const
        {
            return bam && !forceFQ;
        }
        
        FileName index;
        
        // Flip before matching?
        bool flipBefore;
        
        // Flip sequin reads?
        bool flip;
        
        // K-mer length
        KmerLen k;
        
        // How to run Kallisto?
        Product prod;
        
        // How many k-mers to skip?
        Counts skipKM;
        
        // Classification rule (minimum percentage matching)
        Proportion rule;
        
        // Number of threads
        Counts thr;
        
        // Only "sequin" and "ladder" bins?
        bool onlySeqLad;
        
        bool forceFQ;
        
        // BAM mode?
        bool bam;
        
        bool report;
    };
}

#endif
