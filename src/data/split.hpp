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
        SOptions() : skipMerge(false), calibS(-1), calibL(-1) {}
        
        bool skipMerge;
        
        // How much to calibrate for sequins and ladder?
        Proportion calibS, calibL;
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

        // Number of sample reads before and after calibration
        Counts bSam, aSam;
        
        // Number of sequin reads before and after calibation
        Counts bSeq, aSeq;
        
        // Output file names from calibration
        FileName o1, o2;
    };

    // Version number by translation
    std::string SVersion(const Reference &, const KStats &);
    
    /*
     * Anaquin is required to take BAM input files. Unfortuantely, the BAM file may not be sorted. Assume sequin reads
     * are relatively small, we should be able to keep them in memory.
     */
    
    void SCombine(const FileName &, SStats &, const SOptions &,
                  const std::map<ChrID, DIntervals<>> * = nullptr,
                  bool shouldDecoyReverseC = false);
    
    void writeKmers(const FileName &, const SStats &, const SOptions &);

    // Write R for somatic ladder by k-mers
    void writeSomatic(const FileName &, const FileName &, const SOptions &);
    
    void SMerge(SStats &stats, const SOptions &o, const std::set<Bin> &only);

    inline void SWriteLadder(std::shared_ptr<Ladder> l3, const FileName &file, const SStats &stats, const SOptions &o)
    {
        const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
        
        o.generate(file);
        o.writer->open(file);
        o.writer->write((boost::format(format) % "Name" % "Stoch" % "Unit" % "Med" % "Read").str());
        
        for (const auto &seq : stats.K.seqs)
        {
            if (isSubstr(seq, "LD_"))
            {
                auto write = [&](const SequinID &seq)
                {
                    const auto u = stats.R.d2u.meds.count(seq) ? stats.R.d2u.meds.at(seq) : NAN;
                    const auto m = u;
                    const auto x = stats.K.sqc.count(seq) ? stats.K.sqc.at(seq) : 0;
                    
                    if (l3->contains(seq))
                    {
                        o.writer->write((boost::format(format) % seq
                                                               % 1
                                                               % l3->input(seq)
                                                               % (std::isnan(m) ? MISSING : toString(m, 2))
                                                               % x).str());
                    }
                    else
                    {
                        o.writer->write((boost::format(format) % seq
                                                               % MISSING
                                                               % MISSING
                                                               % (std::isnan(m) ? MISSING : toString(m, 2))
                                                               % x).str());
                    }
                };
                
                write(seq);
            }
        }
        
        o.writer->close();
    };

    void SWriteLCopy     (const FileName &, const SOptions &);
    void SWriteLDensity  (const FileName &, const SOptions &o);
    void SWriteLVariation(const FileName &, const SOptions &o);

    template <typename Stats> void SWriteLadderPostCalib(std::shared_ptr<Ladder> l3,
                                                         const Stats &stats,
                                                         const SOptions &o)
    {
        writeKmers(o.name + "_calibrated_kmers.tsv", stats, o);
        
        // Synthetic ladder after calibration
        SWriteLadder(l3, o.name + "_ladder_calibrated.tsv", stats, o);
        
        SWriteLCopy(o.name + "_ladderCopy.R", o);
        SWriteLDensity(o.name + "_ladderDensity.R", o);
    }

    /*
     * Synethetic calibration by calibrating to the lowest abundant sequin. Read pre-calibration from
     * "src" and write calibrated BAM or FASTQ to "o".
     */
    
    template <typename O, typename S, typename F> GSynthetic::Results
            SCalibSynthetic(const S &s,
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
            const auto rr = GSynthetic::calibrate(s.K, src, o, o, o.name);

            auto o2 = o;
            o2.work = tmpPath();
            o2.logger = std::shared_ptr<Writer<>>(new MockWriter());
            o2.output = nullptr; // Suppress messages for the second analysis in f()

            if (o.bam && o.forceFQ)
            {
                o2.bam = false; // The inputs are FQ
            }
            
            // Run analytics on calibrated synthetic (but not writing to the original directory)
            const auto tmp = f(rr.r1, rr.r2, o2); removeD(o2.work);

            o.info("Completed ladder calibration");
            SWriteLadderPostCalib(l3, tmp, o);
            
            return rr;
        }
        catch (const std::runtime_error &ex)
        {
            o.warn("Ladder calibration failed.");
            o.warn(ex.what());
            return GSynthetic::Results();
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
