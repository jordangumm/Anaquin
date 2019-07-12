#include <htslib/bgzf.h>
#include "tools/calibrate.hpp"
#include "Genomics/Genomics.hpp"
#include "writers/fq_writer.hpp"
#include "parsers/parser_fq.hpp"
#include "writers/bam_writer.hpp"
#include "parsers/parser_bam.hpp"

using namespace Anaquin;

struct BAMCalibrator : public Calibrator
{
    virtual ~BAMCalibrator() {}
    
    Result calibrate(const KStats &, Probability p, const WriterOptions &wo) override
    {
        RandomSelection r(1.0 - p);

        BAMWriter w; w.open(o);
        Result rr; rr.o1 = o;
        
        std::set<ReadName> ns;
        
        ParserBAM::parse(f, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
        {
            if (i.p && !(i.p % 10000)) { wo.wait(toString(i.p)); }
            x.lName();
            
            // Trimmed name
            const auto rn = trimRName(x.name);
            
            // Calibrate both pair-ended as they use the same trimmed read name
            if (r.select(rn))
            {
                // Only increment for a paired (not each read)
                if (!ns.count(rn))
                {
                    rr.n++;
                }
                
                ns.insert(rn);
                w.write(x);
            }
        });

        w.close();
        return rr;
    }
    
    Result calibrate(const KStats &stats, Selection &r, const WriterOptions &wo) override
    {
        Result rr; rr.o1 = o;
        BAMWriter w; w.open(o);
        
        ParserBAM::parse(f, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
        {
            if (i.p && !(i.p % 10000)) { wo.wait(toString(i.p)); }
            x.lName();

            const auto r_ = trimRName(x.name);
            assert(stats.r1.at(LD).count(r_));
            
            const auto s1 = GSeq2Std(stats.r1.at(LD).at(r_));
            assert(r.count(s1));
            
            if (r[s1]->select(r_))
            {
                rr.n++;
                w.write(x);
            }
        });
        
        w.close();
        return rr;
    }
    
    FileName f, o;
};

struct FQCalibrator : public Calibrator
{
    virtual ~FQCalibrator() {}

    Result calibrate(const KStats &, Probability p, const WriterOptions &wo) override
    {
        RandomSelection r(1.0 - p);
        
        Result rr;
        rr.o1 = o1; rr.o2 = o2;

        auto w1 = FQWriter(o1);
        auto w2 = FQWriter(o2);
        auto i = 0;
        
        ParserFQ::parse(Reader(f1), Reader(f2), [&](const ParserFQ::Data &x)
        {
            if (i && !(i % 10000)) { wo.wait(toString(i)); } i++;

            auto n1 = trimRName(x.name1);
            auto n2 = trimRName(x.name2);
            assert(n1 == n2);
            
            if (r.select(n1))
            {
                rr.n++;
                w1.write(n1, x.seq1, x.qual1);
                w2.write(n2, x.seq2, x.qual2);
            }
        });
        
        w1.close(); w2.close();
        return rr;
    }
            
    Result calibrate(const KStats &stats, Selection &r, const WriterOptions &wo) override
    {
        Result rr;
        rr.o1 = o1;
        rr.o2 = o2;
        
        auto i = 0;
        std::ofstream w1, w2;
        w1.open(o1); w2.open(o2);

        ParserFQ::parse(Reader(f1), Reader(f2), [&](const ParserFQ::Data &x)
        {
            if (i && !(i % 10000)) { wo.wait(toString(i)); } i++;

            const auto r1 = trimRName(x.name1);
            const auto r2 = trimRName(x.name2);
            
            assert(r1 == r2);
            assert(stats.r1.at(LD).count(r1));
            assert(stats.r2.at(LD).count(r2));
            
            const auto s1 = GSeq2Std(stats.r1.at(LD).at(r1));
            const auto s2 = GSeq2Std(stats.r1.at(LD).at(r2));
            assert(r.count(s1) && r.count(s2));
            
            auto write = [&]()
            {
                w1 << compressGZ("@" + x.name1 + "\n" + x.seq1 + "\n" + x.opt1 + "\n" + x.qual1 + "\n");
                w2 << compressGZ("@" + x.name1 + "\n" + x.seq1 + "\n" + x.opt1 + "\n" + x.qual1 + "\n");
            };
            
            if (s1 != s2 || r[s1]->select(r1))
            {
                rr.n++;
                write();
            }
        });

        w1.close(); w2.close();
        return rr;
    }
            
    FileName f1, f2, o1, o2;
};

std::shared_ptr<Calibrator> Calibrator::createBAM(const FileName &f, const FileName &o)
{
    auto x = std::shared_ptr<BAMCalibrator>(new BAMCalibrator());
    x->f = f; x->o = o;
    return x;
}

std::shared_ptr<Calibrator> Calibrator::createFQ(const FileName &f1, const FileName &f2,
                                                 const FileName &o1, const FileName &o2)
{
    auto x = std::shared_ptr<FQCalibrator>(new FQCalibrator());
    x->f1 = f1; x->f2 = f2; x->o1 = o1; x->o2 = o2;
    return x;
}
