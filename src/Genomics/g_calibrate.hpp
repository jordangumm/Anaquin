#ifndef G_CALIBRATE_HPP
#define G_CALIBRATE_HPP

#include "Genomics/g_split.hpp"
#include "parsers/parser_bambed.hpp"

namespace Anaquin
{
    typedef std::map<ChrID, std::map<Locus, Proportion>> NormFactors;

    struct GCalibrate
    {
        enum class Method
        {
            Mean,
            Median,
            Percent
        };
        
        struct SampledInfo
        {
            SequinID rID;
            
            // Alignment coverage for the endogenous sample
            Coverage endo;
            
            // Alignment coverage before subsampling
            Coverage before;
            
            // Alignment coverage after subsampling
            Coverage after;
            
            // Number of alignments before and after
            Counts nEndo, nBefore, nAfter;
            
            // Normalization factor
            Proportion norm;
        };
        
        struct CalibrateStats
        {
            std::set<ReadName> trimmed;
            
            Counts nh = 0; // Total number of human reads (including calibration regions)
            Counts nd = 0; // Total number of sequin reads (including calibration regions)

            // Number of sequin reads outside calibration region (before calibration)
            Counts nSeqOut = 0;
            
            // Interval statistics for human reads before calibration
            ParserBAMBED::Stats es;
            
            // Interval statistics for sequin reads before calibration
            ParserBAMBED::Stats ss;
            
            std::vector<double> allBeforeEndoC;
            std::vector<double> allBeforeSeqsC;
            
            // Required for summary statistics
            std::vector<double> allNorms;

            // Normalization for each region
            NormFactors norms;

            std::map<ChrID, std::map<Locus, SampledInfo>> c2v;
            
            // Summary statistics for normalization factors
            inline double normMean() const
            {
                return SS::mean(allNorms);
            }

            // Summary statistics for normalization factors
            inline double normSD() const
            {
                return SS::SD(allNorms);
            }
            
            // Average sequence coverage for endogenous before normalization
            inline double meanBEndo() const
            {
                return SS::mean(allBeforeEndoC);
            }
            
            // Average sequence coverage for sequins before normalization
            inline double meanBSeqs() const
            {
                return SS::mean(allBeforeSeqsC);
            }            
        };

        struct GenomeSequins
        {
            Counts nEndo = 0;
            Counts nSeqs = 0;
        };
        
        struct Stats
        {
            // Before calibration
            GSplit::Stats S1;
            
            // After calibration
            GSplit::Stats S2;

            CalibrateStats cStats;
            
            Coverage afterSeqs;
            
            GenomeSequins tBefore, tAfter;
            GenomeSequins sBefore, sAfter;
            
            std::map<ChrID, std::map<Locus, SampledInfo>> c2v;
        };
        
        enum class CalibrateMode
        {
            TwoBAM,
            Combined
        };

        struct Options : public GSplit::Options
        {
            Options() : GSplit::Options(), trim(true)
            {
                onlySeqLad = true;
            }

            // How to calculate coverage?
            Method meth = Method::Mean;

            CalibrateMode cMode;
            
            Base edge;

            // Defined only if Method::Percent
            Proportion p = NAN;
            
            // Defined only if meth==Reads
            Counts reads;
            
            // Should we trim alignments before calibration?
            bool trim;
        };

        static GenomeSequins tBefore(const CalibrateStats &);
        static GenomeSequins sBefore(const CalibrateStats &);
        static GenomeSequins tAfter (const CalibrateStats &, const ParserBAMBED::Stats &);
        static GenomeSequins sAfter (const CalibrateStats &, const ParserBAMBED::Stats &);
        
        static double afterSeqsC(const Chr2DInters &,
                                 std::map<ChrID, std::map<Locus, SampledInfo>> &,
                                 const Options &);
        
        static ParserBAMBED::Stats calibrate(const FileName    &,
                                             const FileName    &,
                                             const NormFactors &,
                                             const Chr2DInters &,
                                             const Chr2DInters &,
                                             const std::set<ReadName> &,
                                             CalibrateStats &,
                                             const Options &);
        
        static CalibrateStats twoBAM(const FileName &,
                                     const FileName &,
                                     const Chr2DInters &,
                                     const Chr2DInters &,
                                     const Options &);

        static CalibrateStats combined(const FileName &,
                                       const FileName &,
                                       const Chr2DInters &,
                                       const Chr2DInters &,
                                       const Chr2DInters &,
                                       SStats &,
                                       const Options &);

        static Stats analyze(const FileName &, const FileName &, const Options &);

        static void report(const FileName &, const Options &o = Options());
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };    
}

#endif
