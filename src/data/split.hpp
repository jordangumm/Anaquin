#ifndef SPLIT_HPP
#define SPLIT_HPP

#include <fstream>
#include "Kallisto.hpp"
#include "data/fastq.hpp"
#include "tools/tools.hpp"
#include "stats/stats.hpp"
#include "stats/linear.hpp"
#include "stats/ss/stats.hpp"
#include "stats/analyzer.hpp"
#include "data/resources.hpp"
#include "parsers/parser_bam.hpp"
#include "Genomics/g_synthetic.hpp"

namespace Anaquin
{
    struct SOptions : public KOptions, public AnalyzerOptions<HumanAssembly>
    {
        SOptions() : skipMerge(true), oneS(-1), oneL(-1), secL(-1) {}
        
        bool skipMerge;
        
        // How much to calibrate for sequins and ladder in stage one?
        Proportion oneS, oneL;

        // Stage two ladder calibration
        Proportion secL;

        // Are we doing ladder calibration?
        inline bool isSCalib() const
        {
            return oneS != NO_CALIBRATION;
        }

        // Are we doing ladder calibration?
        inline bool isLCalib() const
        {
            return (oneL != NO_CALIBRATION) || (secL != NO_CALIBRATION);
        }

        // Ladder calibration can only be done once, which one?
        inline Proportion ladP() const
        {
            return (oneL != NO_CALIBRATION ? oneL : secL);
        }
    };
    
    struct SStats
    {
        // Kallisto statistics
        KStats K;
        
        // Dilution
        inline Proportion dil() const { return 1.0 - K.binP(ES); }

        typedef std::map<SequinID, std::vector<double>> KMCount;
        
        struct Abundance
        {
            // Shared k-mer counts
            KMCount s2s;
            
            // Unique k-mer counts
            KMCount s2u;
            
            struct Descriptive
            {
                // Descriptive statistics (unique k-mers)
                std::map<SequinID, Counts> mins, maxs;
                
                // Descriptive statistics (unique k-mers)
                std::map<SequinID, KMCoverage> mus, q25, q75, meds, sds;
            };
            
            // Descriptive statistics for shared and unique k-mers
            Descriptive d2s, d2u;
        };
        
        Abundance R;

        // Number of human reads for combined (hacked for combined)
        Counts hackHumanReads = 0;
    };
    
    // Statistics for calibration
    struct SCStats
    {
        Proportion p;
        
        // Target number of reads
        Counts tar;

        Counts bSam = 0; // Sample reads before calibration
        Counts aSam = 0; // Sample reads after calibration        
        Counts bSeq = 0; // Sequin reads before calibration
        Counts aSeq = 0; // Sequin reads after calibration
        
        // Output file names from calibration
        FileName o1, o2;
    };

    // Version number by translation
    std::string SVersion(const Reference &, const KStats &);
    
    struct SAlignStats
    {
        SAlignStats()
        {
            ins[Bin::ES]; // Make sure we always have an entry
            ins[Bin::GR]; // Make sure we always have an entry
        }

        // Average insertion size
        double mInsert(Bin b = Bin::GR, Base min = 10, Base max = 1000) const;
        
        // Frequency table for insertion size of each bin
        std::map<Bin, std::map<Base, Counts>> ins;
    };
    
    /*
     * Anaquin is required to take BAM input files. Unfortuantely, the BAM file may not be sorted.
     * Assume sequin reads are relatively small, we should be able to keep them in memory.
     */
    
    SAlignStats SAlign(const FileName &, SStats &, const SOptions &,
                       const std::map<ChrID, DIntervals<>> * = nullptr,
                       bool shouldDecoyReverseC = false);
    
    void writeLTable(const FileName &, const FileName &, const SOptions &);
    void writeSTable(const FileName &, const FileName &, const SOptions &, Counts nExp, Counts nObs, Counts nCV, const Label &, const Label &);
    
    void writeInsertR(const FileName &, const FileName &, const SOptions &);
    void writeInsertTSV(const FileName &, const SAlignStats &, const SOptions &);
    
    void writeKmers(const FileName &, const SStats &, const SOptions &);

    // Write R for somatic ladder by k-mers
    void writeSomatic(const FileName &, const FileName &, const SOptions &);
    
    void SMerge(SStats &stats, const SOptions &o, const std::set<Bin> &only);

    void SWriteLadder(std::shared_ptr<Ladder>, const FileName &, const SStats &, const SOptions &);

    void SWriteLCopy(const FileName &, const Label &, const SOptions &);
    void SWriteLDensity(const FileName &, const Label &, const Label &, const SOptions &o);

    template <typename Stats> void SWriteLadderPostCalib(std::shared_ptr<Ladder> l3,
                                                         const Stats &stats,
                                                         const SOptions &o)
    {
        const auto l1 = (!o.isLCalib() ? "_kmers.tsv"  : "_kmers_calibrated.tsv");
        const auto l2 = (!o.isLCalib() ? "_ladder.tsv" : "_ladder_calibrated.tsv");

        // Generating "_kmers.tsv" or "_kmers_calibrated.tsv"
        writeKmers(o.name + l1, stats, o);
        
        // Generating "_ladder.tsv" or "_ladder_calibrated.tsv"
        SWriteLadder(l3, o.name + l2, stats, o);
        
        if (o.report)
        {
            // R-script for ladder abundance
            SWriteLCopy("report_files/" + o.name + "_ladderCopy.R", l2, o);
            
            // R-script for ladder density
            SWriteLDensity("report_files/" + o.name + "_ladderDensity.R", l1, l2, o);
        }
    }

    /*
     * Synethetic calibration by calibrating to the lowest abundant sequin. Read pre-calibration from
     * "src" and write calibrated BAM or FASTQ to "o".
     */
    
    template <typename O, typename S, typename F> LadderInternal::Results
            StageTwoLadderCalibration(const S &s,
                            F f,
                            const Path &src,
                            const O &o,
                            std::shared_ptr<Ladder> l3)
    {
        assert(l3);
        
        // Generating "_ladder.tsv" before calibration
        SWriteLadder(l3, o.name + "_ladder.tsv", s, o);
        
        try
        {
            // Apply stage two internal ladder calibration
            const auto rr = LadderInternal::calibrate(s.K, src, o, o.secL, o, o.name);

            auto o2 = o;
            o2.work = tmpPath();
            o2.logger = std::shared_ptr<Writer<>>(new MockWriter());
            o2.output = nullptr; // Suppress messages for the second analysis in f()

            if (o.bam && o.forceFQ)
            {
                o2.bam = false; // The inputs are FQ
            }
            
            // Run analytics on calibrated ladders (but not writing to the original directory)
            const auto tmp = f(rr.r1, rr.r2, o2); removeD(o2.work);

            SWriteLadderPostCalib(l3, tmp, o);
            
            return rr;
        }
        catch (const std::runtime_error &ex)
        {
            o.warn("Ladder calibration failed.");
            o.warn(ex.what());
            return LadderInternal::Results();
        }
    }

    // Writing reads to TSV file
    void SWriteReads(Product, const FileName &, const SStats &, const SOptions &);

    void SKallisto(SStats &, const FileName &, const FileName &, const SOptions &);

    /*
     * -------------------- Calibration by percentage --------------------
     */

    struct SCalibratePFiles
    {
        virtual FileName i1(const SOptions &) const = 0;
        virtual FileName i2(const SOptions &) const = 0;
        virtual FileName o1(const SOptions &) const = 0;
        virtual FileName o2(const SOptions &) const = 0;
    };
    
    /*
     * Default implementation for SCalibratePFiles
     */
    
    struct SCalibrateDefault : public SCalibratePFiles
    {
        SCalibrateDefault(const Label &x, const Path &f1 = "", const Path &f2 = "") : x(x), f1(f1), f2(f2) {}
        
        FileName i1(const SOptions &o) const override
        {
            const auto f1_ = f1.empty() ? o.work : f1;
            return o.writeBAM() ? f1_ + "/" + o.name + "_" + x + ".bam" :
                                  f1_ + "/" + o.name + "_" + x + "_1.fq.gz";
        }
        
        FileName i2(const SOptions &o) const override
        {
            const auto f1_ = f1.empty() ? o.work : f1;
            return o.writeBAM() ? "" : f1_ + "/" + o.name + "_" + x + "_2.fq.gz";
        }
        
        FileName o1(const SOptions &o) const override
        {
            const auto f2_ = f2.empty() ? o.work : f2;
            return o.writeBAM() ? f2_ + "/" + o.name + "_" + x + "_calibrated.bam" :
                                  f2_ + "/" + o.name + "_" + x + "_calibrated_1.fq.gz";
        }
        
        FileName o2(const SOptions &o) const override
        {
            const auto f2_ = f2.empty() ? o.work : f2;
            return o.writeBAM() ? "" : f2_ + "/" + o.name + "_" + x + "_calibrated_2.fq.gz";
        }
        
        const Label x;
        const FileName f1, f2;
    };
    
    // Calibrate by percentage for a bin. x must be "GR" or "LD".
    SCStats SCalibrateP(Bin, Proportion, const SStats &, const SOptions &, const SCalibratePFiles &);
}

#endif
