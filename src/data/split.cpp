#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include "data/split.hpp"
#include "tools/tools.hpp"
#include "data/reader.hpp"
#include "tools/samtools.hpp"
#include "tools/calibrate.hpp"
#include "writers/r_writer.hpp"
#include "writers/fq_writer.hpp"
#include "writers/bam_writer.hpp"
#include "writers/sam_writer.hpp"
#include "parsers/parser_bam.hpp"

extern std::vector<const Anaquin::ReadName *> __HackBAMDecoy2__;
extern std::map<std::string, std::map<int, Anaquin::KMInfo>> __KMInfo__;

// Defined in ProcessReads
extern Anaquin::AlignedReads __CombinedBAM1__;

// Defined in ProcessReads
extern bam_hdr_t *__CombinedBAMHeader__;

using namespace Anaquin;

// Defined in Kallisto
extern KStats Kallisto(const FileName &, const FileName &, const KOptions &);

template <typename T1, typename T2> void SKFreq(const T1 &x, T2 &m)
{
    for (const auto &i : x) // For all sequins
    {
        for (const auto &j : i.second) // For all k-mers in the sequin
        {
            assert(j.second);
            m[i.first].push_back(j.second);
        }
    }
};

template <typename T> void SQuantKM(T &x)
{
    // Descriptive statistics for unique k-mers
    for (auto &i : x.s2u)
    {
        std::sort(i.second.begin(), i.second.end());
        x.d2u.mus[i.first]  = SS::mean(i.second);
        x.d2u.q25[i.first]  = SS::hackQuant(i.second, 0.25);
        x.d2u.q75[i.first]  = SS::hackQuant(i.second, 0.75);
        x.d2u.sds[i.first]  = SS::SD(i.second);
        x.d2u.mins[i.first] = *(i.second.begin());
        x.d2u.meds[i.first] = SS::hackQuant(i.second, 0.50);
        x.d2u.maxs[i.first] = *(i.second.rbegin());
    }
    
    // Descriptive statistics for shared k-mers
    for (auto &i : x.s2s)
    {
        std::sort(i.second.begin(), i.second.end());
        x.d2u.mus[i.first]  = SS::mean(i.second);
        x.d2s.q25[i.first]  = SS::hackQuant(i.second, 0.25);
        x.d2s.q75[i.first]  = SS::hackQuant(i.second, 0.75);
        x.d2s.sds[i.first]  = SS::SD(i.second);
        x.d2s.mins[i.first] = *(i.second.begin());
        x.d2s.meds[i.first] = SS::hackQuant(i.second, 0.50);
        x.d2s.maxs[i.first] = *(i.second.rbegin());
    }
};

void Anaquin::SWriteReads(Product mode, const FileName &file, const SStats &stats, const SOptions &o)
{
#ifdef WRITE_READS
    const auto format = "%1%\t%2%\t%3%\t%4%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Name" % "Read1" % "Read2" % "Label").str());
    
    auto write = [&](const std::map<ReadName, SequinID> &x, const std::map<ReadName, SequinID> &y, const Label &l)
    {
        for (const auto &i : x)
        {
            o.writer->write((boost::format(format) % i.first % i.second % y.at(i.first) % l).str());
        }
    };
    
    forBin([&](Bin c)
    {
        std::string s;
        
        switch (mode)
        {
            case Product::Genomics:
            {
                switch (c)
                {
                    case HP: { s = "HP";         break; }
                    case IF: { s = "Info";       break; }
                    case MT: { s = "Mito";       break; }
                    case HL: { s = "HLA";        break; }
                    case LD: { s = "Ladder";     break; }
                    case IM: { s = "Immune";     break; }
                    case ES: { s = "Sample";     break; }
                    case VC: { s = "Vector";     break; }
                    case SV: { s = "Structural"; break; }
                    case GR: { s = "Germline";   break; }
                    case SO: { s = "Somatic";    break; }
                    case MS: { s = "Micro";      break; }
                    case MI: { s = "MSI";        break; }
                }
                
                break;
            }

            case Product::Meta:
            {
                switch (c)
                {
                    case IF: { s = "Info";    break; }
                    case LD: { s = "Ladder";  break; }
                    case ES: { s = "Sample";  break; }
                    case VC: { s = "Vector";  break; }
                    case GR: { s = "Sequins"; break; }
                    default: { break; }
                }

                break;
            }

            case Product::RNA:
            {
                switch (c)
                {
                    case IF: { s = "Info";    break; }
                    case ES: { s = "Sample";  break; }
                    case GR: { s = "Sequins"; break; }
                    default: { break; }
                }
                
                break;
            }

            default: { break; }
        }

        if (stats.K.r1.count(c))
        {
            write(stats.K.r1.at(c), stats.K.r2.at(c), s);
        }
    });

    o.writer->close();
#endif
}

/*
 * Post-calculation from completed Kallisto statistics. This function should be called after Kallisto. It
 * is responsible for constructing various ladders.
 */

static void SPostCal(SStats &stats, const SOptions &o)
{
    o.info("Quantifying k-mers");

    SKFreq(stats.K.uniqs, stats.R.s2u);
    SKFreq(stats.K.shared, stats.R.s2s);
    SQuantKM(stats.R);
     
    if (stats.R.d2s.mins.empty() && stats.R.d2u.mins.empty())
    {
        o.warn("Sequin not found. Please check your input files.");
    }
    
    o.info("Merging k-mers");
    std::set<Bin> only;
    
    switch (o.prod)
    {
        case Product::RNA:  { only = std::set<Bin> { ES, GR, VC };         break; }
        case Product::Meta: { only = std::set<Bin> { ES, IF, GR, LD, VC }; break; }
        case Product::Genomics:
        {
            only = std::set<Bin> { IF, MT, MS, HL, HP, LD, SV, IM, ES, GR, SO, MI, VC };
            break;
        }
    }

    /*
     * Never generate sample BAM if the input source is BAM as they could be quite large.
     */
    
    if (o.bam)
    {
        assert(only.count(ES));
        only.erase(only.find(ES));
    }

    // Merge and generate FASTQ for each bin
    SMerge(stats, o, only);

    assert(stats.dil() >= 0 && stats.dil() <= 1.0);
}

struct AbstractMerge
{
    virtual void init(const SOptions &o, Bin m, const Label &x) = 0;
    virtual void close() = 0;
    virtual void merge(SStats &, Bin, Thread, const SOptions &) = 0;
};

struct FQMerge : public AbstractMerge
{
    void init(const SOptions &o, Bin m, const Label &x) override
    {
        const auto f1_ = o.work + "/" + o.name + "_" + x + "_1.fq.gz";
        const auto f2_ = o.work + "/" + o.name + "_" + x + "_2.fq.gz";

        auto i1 = std::find_if(f1.begin(), f1.end(), [&](const std::pair<Bin, FileName> &i)
        {
            return i.second == f1_;
        });

        auto i2 = std::find_if(f2.begin(), f2.end(), [&](const std::pair<Bin, FileName> &i)
        {
            return i.second == f2_;
        });

        f1[m] = f1_;
        f2[m] = f2_;

        if (i1 != f1.end())
        {
            assert(i2 != f2.end());
            o1[m] = o1[i1->first];
            o2[m] = o2[i2->first];
            prim[m] = false;
        }
        else
        {
            prim[m] = true;
            o1[m] = std::shared_ptr<std::ofstream>(new std::ofstream(f1[m], std::ios::binary | std::ios::out));
            o2[m] = std::shared_ptr<std::ofstream>(new std::ofstream(f2[m], std::ios::binary | std::ios::out));
        }
    }
    
    void close() override
    {
        for (auto &i : o1) { if (prim[i.first]) { i.second->close(); } }
        for (auto &i : o2) { if (prim[i.first]) { i.second->close(); } }
    }

    void merge(SStats &stats, Bin b, Thread i, const SOptions &o) override
    {
        assert(f1.count(b) && f2.count(b));
        assert(stats.K.f1[b][i].size() == stats.K.f2[b][i].size());
        
        auto myCopyGZ = [&](const FileName &src, std::shared_ptr<std::ofstream> w)
        {
            switch (copyGZ(src, w))
            {
                case CopyGZStatus::Failed:
                {
                    o.logWarn(src + " malformed");
                    return false;
                }
                    
                case CopyGZStatus::Corrected:
                {
                    o.logWarn(src + " malformed but corrected");
                    return true;
                }

                case CopyGZStatus::Success: { return true; }
            }
        };

        if (!myCopyGZ(stats.K.f1[b][i], o1[b])) { return; }
        if (!myCopyGZ(stats.K.f2[b][i], o2[b])) { return; }

        removeF(stats.K.f1[b][i]); // Remove the old file
        removeF(stats.K.f2[b][i]); // Remove the old file
    }
    
    std::map<Bin, bool> prim;
    std::map<Bin, FileName> f1, f2;
    std::map<Bin, std::shared_ptr<std::ofstream>> o1, o2;
};

struct BAMMerge : public AbstractMerge
{
    void init(const SOptions &o, Bin m, const Label &x) override
    {
        // BAM file writing to
        const auto file = o.work + "/" + o.name + "_" + x + ".bam";
        
        auto iter = std::find_if(f.begin(), f.end(), [&](const std::pair<Bin, FileName> &i)
        {
            return i.second == file;
        });

        // Reuse the same pointer?
        if (iter != f.end())
        {
            w[m] = w[iter->first];
            prim[m] = false;
        }
        else
        {
            prim[m] = true;
            assert(__CombinedBAMHeader__);
            w[m] = bgzf_open(file.c_str(), "w");
            bam_hdr_write(w[m], __CombinedBAMHeader__);
        }

        f[m] = file;
    }
    
    void close() override
    {
        for (auto &i : w) { if (prim[i.first]) { bgzf_close(i.second); } }
    }
    
    void copyBAM(const FileName &src, BGZF *dst)
    {
        ParserBAM::parse(src, [&](const ParserBAM::Data &x, const ParserBAM::Info &)
        {
            bam_write1(dst, (const bam1_t *) x.b());
        });
    }

    void merge(SStats &stats, Bin b, Thread i, const SOptions &) override
    {
        assert(f.count(b) && w.count(b));
        copyBAM(stats.K.f1[b][i], w[b]);
        removeF(stats.K.f1[b][i]); // Remove the old file
    }

    std::map<Bin, bool> prim;
    std::map<Bin, FileName> f;
    std::map<Bin, BGZF *> w;
};

void Anaquin::SMerge(SStats &stats, const SOptions &o, const std::set<Bin> &only)
{
    auto xs = std::vector<std::shared_ptr<AbstractMerge>>();
    
    if (o.writeBAM()) { xs.push_back(std::shared_ptr<AbstractMerge>(new BAMMerge())); }
    else              { xs.push_back(std::shared_ptr<AbstractMerge>(new FQMerge()));  }
    
    auto bin2Str = [&](Bin x)
    {
        switch (x)
        {
            case SV: { return "sv";     }
            case IF: { return "info";   }
            case HL: { return "hla";    }
            case MT: { return "mito";   }
            case LD: { return "ladder"; }
            case IM: { return "immune"; }
            case ES: { return "sample"; }
            case VC: { return "vector"; }
            case HP:
            case MS:
            case GR:
            case MI:
            case SO: { return "sequin"; }
        }
    };
    
    forBin([&](Bin x)
    {
        if (!only.count(x))
        {
            return;
        }
        else if (o.onlySeqLad && x != GR && x != LD)
        {
            return;
        }
        
        for (auto &i : xs)
        {
            i->init(o, x, bin2Str(x));
        }
    });
    
    // For each category...
    for (auto c : stats.K.f1)
    {
        auto bin = c.first;

        if (!only.count(bin))
        {
            continue;
        }
        else if (o.onlySeqLad && bin != GR && bin != LD)
        {
//            continue;
        }

        auto m = std::map<Bin, std::string>
        {
            { ES, "ES" },
            { IF, "IF" },
            { MT, "MT" },
            { MS, "MS" },
            { HL, "HL" },
            { HP, "HP" },
            { LD, "LD" },
            { SV, "SV" },
            { IM, "IM" },
            { ES, "ES" },
            { GR, "GR" },
            { SO, "SO" },
            { VC, "VC" },
        };
        
        // For each thread...
        for (auto i = 0u; i < stats.K.f1[bin].size(); i++)
        {
            o.info("Merging thread " + toString(i));

            for (auto &j : xs)
            {
                j->merge(stats, bin, i, o);
            }
        }
    }
    
    for (auto &i : xs) { i->close(); }
    
    assert(!stats.K.work.empty());
    
    // Remove the directory keeping Kallisto files
    removeD(stats.K.work);
}

SCStats Anaquin::SCalibrateP(Bin b, Proportion p, const SStats &stats, const SOptions &o, const SCalibratePFiles &f)
{
    assert(b == GR || b == LD);
    o.info("Calibrating by percentage");

    SCStats x;
    assert(o.prod == Product::Meta || o.prod == Product::RNA);
    
    // Number of sample paired reads (equal before and after calibration)
    x.bSam = x.aSam = stats.K.binN(ES);

    // That's what we need to exclude
    b = (b == LD ? GR : LD);
    
    // Number of paired reads before calibration
    x.bSeq = stats.K.total() - x.bSam - stats.K.binN(b) - stats.K.binN(VC) - stats.K.binN(IF);

    if (x.bSam == 0) { o.warn("Number of sample reads is zero. Scaling factor set to 1."); }
    if (x.bSeq == 0) { o.warn("Number of sequin reads is zero. Scaling factor set to 0."); }

    if (p == -1)
    {
        o.info("Skipped calibration");
        x.p  = -1;
        x.o1 = f.i1(o); // No calibration
        x.o2 = f.i2(o); // No calibration
        x.aSeq = x.bSeq;
        return x;
    }
    
    // Number of target sequin reads after calibration
    x.tar = (p / (1.0 - p)) * x.bSam;
    
    // Make sure our target doesn't goto zero if the sample is non-zero
    if (x.bSam && !x.tar)
    {
        x.tar = x.bSam;
    }
    
    // Percentage for calibration
    x.p = (x.bSam == 0) ? 1.0 : (x.bSeq == 0) ? 0.0 : (x.tar >= x.bSeq ? 1.0 : ((float) x.tar) / x.bSeq);
    
    assert(x.p >= 0.0 && x.p <= 1.0);   
    
    auto c = (o.writeBAM()) ? Calibrator::createBAM(f.i1(o), f.o1(o)) : Calibrator::createFQ(f.i1(o), f.i2(o), f.o1(o), f.o2(o));
    auto r = c->calibrate(stats.K, x.p, o);
    
    x.o1   = r.o1;
    x.o2   = r.o2;
    x.aSeq = r.n;
    assert(x.bSeq >= x.aSeq);

    return x;
}

void Anaquin::SWriteLCopy(const FileName &file, const SOptions &o)
{
    const auto p1 = o.calibL ? "_ladder_calibrated.tsv" : "_ladder.tsv";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotKSynthetic()) % date()
                                                     % o.cmd
                                                     % o.work
                                                     % (o.name + p1)).str());
    o.writer->close();
}

void Anaquin::SWriteLDensity(const FileName &file, const SOptions &o)
{
    const auto p1 = "_calibrated_kmers.tsv";
    const auto p2 = "_ladder_calibrated.tsv";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotLDensity()) % date()
                                                   % o.cmd
                                                   % o.work
                                                   % (o.name + p1)
                                                   % (o.name + p2)).str());
    o.writer->close();
}

void Anaquin::SWriteLVariation(const FileName &file, const SOptions &o)
{
    const auto p1 = "_calibrated_kmers.tsv";
    const auto p2 = "_ladder_calibrated.tsv";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotLVariation()) % date()
                                                     % o.cmd
                                                     % o.work
                                                     % (o.name + p1)
                                                     % (o.name + p2)).str());
    o.writer->close();
}

void Anaquin::SKallisto(SStats &stats, const FileName &f1, const FileName &f2, const SOptions &o)
{
    // Running Kallisto
    stats.K = Kallisto(f1, f2, o);
    
    // Post-calculation
    SPostCal(stats, o);

    if (o.bam && stats.hackHumanReads)
    {
        // Sample reads are forced to zero
        assert(stats.K.c1[ES] == 0 && stats.K.c2[ES] == 0);
        
        // Force the human reads
        stats.K.c1[ES] = stats.K.c2[ES] = stats.hackHumanReads;
        
        assert(!std::isnan(stats.K.c1[ES]) && !std::isnan(stats.K.c2[ES]));
    }
}

void Anaquin::writeSomatic(const FileName &file, const FileName &src, const SOptions &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write(RWriter::createLinear(src, o.work, "Allele Frequency Ladder", "Expected Allele Frequency (log2)", "Measured Allele Frequency (log2)", "data$ExpFreq", "data$ObsFreq"));
    o.writer->close();
}

std::string Anaquin::SVersion(const Reference &r, const KStats &x)
{
    std::map<SequinID, Counts> c;
    
    for (const auto &i : x.sqc)
    {
        if (r.t1()->count(i.first))
        {
            c[i.first] = i.second;
        }
    }
    
    return !c.empty() ? r.t1()->at(max(c)) : "-";
}

void Anaquin::writeKmers(const FileName &file, const SStats &stats, const SOptions &o)
{
    /*
     * Consturct a convenient structure for mapping
     */
    
    std::map<std::string, std::map<Kmer, KMInfo *>> m1;

    for (auto &i : __KMInfo__)
    {
        for (auto &j : i.second)
        {
            m1[i.first][j.second.kmer] = &j.second;
        }
    }

    auto add = [&](const SequinKM &x)
    {
        for (const auto &i : x)
        {
            for (const auto &j : i.second)
            {
                const auto k = m1[i.first].count(j.first) ? j.first : revcomp(j.first);
                assert(m1[i.first][k]);
                m1[i.first][k]->abund += j.second;
            }
        }
    };
    
    add(stats.K.uniqs);
    add(stats.K.shared);
    
    const auto format = "%1%\t%2%\t%3%\t%4%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "Sequin"
                                           % "Sequence"
                                           % "Position"
                                           % "Count").str());

    for (const auto &i : __KMInfo__)
    {
        const auto &seq = i.first;
        
        for (const auto &j : i.second)
        {
            o.writer->write((boost::format(format) % seq
                                                   % j.second.kmer
                                                   % j.first
                                                   % j.second.abund).str());
        }
    }

    o.writer->close();
}

void Anaquin::SCombine(const FileName &file, SStats &stats, const SOptions &o, const std::map<ChrID, DIntervals<>> *r1, bool shouldDecoyReverseC)
{
    if (o.bam)
    {
        o.analyze(file);
        __CombinedBAM1__.clear(); __HackBAMDecoy2__.clear();
        
        BAMWriter w;
        
        if (!o.onlySeqLad)
        {
            w.open(o.work + "/" + o.name + "_sample.bam");
        }

        ParserBAM::parse_(file, __CombinedBAM1__, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
        {
            ParserBAM::ParseResult rr;
            
            if (i.p && !(i.p % 10000)) { o.wait(toString(i.p)); }
         
            if (i.p == 0)
            {
                // Required when constructing the header later...
                __CombinedBAMHeader__ = (bam_hdr_t *) x.copyH();
                
                if (!o.onlySeqLad)
                {
                    // Write headers in case there's no sample read
                    w.writeH(x);
                }
            }
            
            // Decoy chromosome?
            const auto decoy = x.mapped && isDecoy(x.cID);
            
            // Fall into sequin region? (only if provided)
            const auto isReg = r1 && r1->count(x.cID) && r1->at(x.cID).overlap(x.l);

            // This read will be selected to Kallisto
            rr.isSeq = decoy || isReg;
            
            if (rr.isSeq)
            {
                x.lSeq(); x.lQual(); x.lName();
            }

            if (shouldDecoyReverseC && decoy && x.isReverseC)
            {
                rr.rc = true;
            }

            /*
             * We don't assume sorted BAM, and there is no way to determine if both paired-end reads
             * are present. We will only count on the first pair.
             */
            
            if (!decoy && x.isFirstPair)
            {
                stats.hackHumanReads++;
            }

            if (!isDecoy(x.cID) && !o.onlySeqLad)
            {
                w.write(x);
            }
            
            return rr;
        });

        w.close();        
        __HackBAMDecoy2__.clear();
        for (const auto &x : __CombinedBAM1__) { __HackBAMDecoy2__.push_back(&(x.first)); }
    }
}
