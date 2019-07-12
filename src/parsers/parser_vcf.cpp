#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "tools/tools.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

static bool isStrelka(const Variant &x)
{
    auto isSNP = [&]()
    {
        return x.fi.count("AU_2_1") &&
               x.fi.count("CU_2_1") &&
               x.fi.count("GU_2_1") &&
               x.fi.count("TU_2_1") &&
               x.fi.count("AU_2_2") &&
               x.fi.count("CU_2_2") &&
               x.fi.count("GU_2_2") &&
               x.fi.count("TU_2_2");
    };
    
    auto isInd = [&]()
    {
        return x.fi.count("TAR_1_1") &&
               x.fi.count("TAR_1_2") &&
               x.fi.count("TAR_2_1") &&
               x.fi.count("TAR_2_2") &&
               x.fi.count("TIR_1_1") &&
               x.fi.count("TIR_1_2") &&
               x.fi.count("TIR_2_1") &&
               x.fi.count("TIR_2_2");
    };
    
    return isSNP() || isInd();
}

std::map<ChrID, Base> ParserVCF::contigs(const Reader &r)
{
    std::string line;
    std::map<ChrID, Base> x;
  
    while (r.nextLine(line))
    {
        if (line[0] != '#')
        {
            break;
        }
        else if (isSubstr(line, "contig"))
        {
            /*
             * Eg: contig=<ID=chr12,length=133275309>
             */
            
            std::vector<std::string> t1, t2;
            
            split(line, "=", t1);
            split(t1[2], ",", t2);

            const auto cID = CHROM(t2[0]);
            split(std::string(t1[3]), ">", t1);

            // 133275309
            x[cID] = stoi(t1[0]);
        }
    }

    return x;
}

void ParserVCF::parse(const Reader &r, Functor f)
{
    htsFile *fp = bcf_open(r.src().c_str(), "r");
    
    if (!fp)
    {
        throw std::runtime_error("Failed to open: " + r.src());
    }
    
    bcf_hdr_t *hdr = bcf_hdr_read(fp);

    if (!hdr)
    {
        throw std::runtime_error("Failed to open: " + r.src());
    }

    int c;

    char  *ps = (char  *) malloc(1024 * sizeof(int));
    int   *pi = (int *)   malloc(1024 * sizeof(int));
    float *pf = (float *) malloc(1024 * sizeof(float));
    
    bcf1_t *line = bcf_init();
    
    while (bcf_read(fp, hdr, line) == 0)
    {
        if (line->n_sample > 2)
        {
            throw std::runtime_error("Number of samples in VCF must not be more than 2.");
        }
        
        bcf_unpack(line, BCF_UN_ALL);
        
        auto run = [&](int al)
        {
            Variant x;

            x.cID  = CHROM(std::string(bcf_seqname(hdr, line)));
            x.name = line->d.id;
            
            x.l.start = x.l.end = line->pos + 1;
            
            x.ref = std::string(line->d.allele[0]);
            x.alt = std::string(line->d.allele[al]);
            
            if (x.ref.size() != 1)
            {
                x.l.end = x.l.start + x.ref.size() - 1;
            }
            
            if (x.alt == "." || x.alt == "")
            {
                return;
            }
            
            auto infos = [&](const std::string &k)
            {
                if (bcf_get_info_string(hdr, line, k.c_str(), &ps, &c) > 0)
                {
                    x.ifs[k] = ps;
                }
            };
            
            auto infof = [&](const std::string &k)
            {
                if (bcf_get_info_float(hdr, line, k.c_str(), &pf, &c) > 0)
                {
                    x.iff[k] = *pf;
                }
            };
            
            auto infoi = [&](const std::string &k)
            {
                if (bcf_get_info_int32(hdr, line, k.c_str(), &pi, &c) > 0)
                {
                    x.ifi[k] = *pi;
                }
            };
            
            infos("CS");
            infos("RT");
            infos("GC");
            infos("GN");
            infos("GT");
            infof("AF");
            infoi("CP");
            infoi("DP");
            infoi("SVLEN");
            
            infoi("QSI");        // Strelka
            infoi("QSS");        // Strelka
            infof("SomaticEVS"); // Strelka
            
            infof("TLOD"); // MuTect
            
            if (x.iff.count("AF")) { x.allF  = x.iff.at("AF"); }
            
            int32_t *g1 = NULL, g2 = 0;
            const auto gt = bcf_get_genotypes(hdr, line, &g1, &g2);

            if (gt >= 2)
            {
                x.gt = g1[0] == g1[1] ? Genotype::Homozygous : Genotype::Heterzygous;
                free(g1);
            }
            else if (gt == 1)
            {
                x.gt = g1[0] == 0 ? Genotype::Homozygous : Genotype::Heterzygous;
                free(g1);
            }

            auto fi = [&](const std::string &key, const std::string &to, int i)
            {
                if (bcf_get_format_int32(hdr, line, key.c_str(), &pi, &c) > i)
                {
                    x.fi[to] = *(pi+i);
                }
            };
            
            auto ff = [&](const std::string &key, const std::string &to, int i)
            {
                if (bcf_get_format_float(hdr, line, key.c_str(), &pf, &c) > i)
                {
                    x.ff[to] = *(pf+i);
                }
            };
            
            fi("TAR", "TAR_1_1", 0);
            fi("TAR", "TAR_1_2", 1);
            fi("TAR", "TAR_2_1", 2);
            fi("TAR", "TAR_2_2", 3);
            
            fi("TIR", "TIR_1_1", 0);
            fi("TIR", "TIR_1_2", 1);
            fi("TIR", "TIR_2_1", 2);
            fi("TIR", "TIR_2_2", 3);
            
            fi("DP", "DP_1", 0);
            fi("DP", "DP_2", 1);
            ff("AF", "AF_1", 0);
            ff("AF", "AF_2", 1);
            
            fi("AU", "AU_1_1", 0);
            fi("AU", "AU_1_2", 1);
            fi("AU", "AU_2_1", 2); // Tumor Tier-1
            fi("AU", "AU_2_2", 3);
            fi("CU", "CU_1_1", 0);
            fi("CU", "CU_1_2", 1);
            fi("CU", "CU_2_1", 2); // Tumor Tier-1
            fi("CU", "CU_2_2", 3);
            fi("GU", "GU_1_1", 0);
            fi("GU", "GU_1_2", 1);
            fi("GU", "GU_2_1", 2); // Tumor Tier-1
            fi("GU", "GU_2_2", 3);
            fi("TU", "TU_1_1", 0);
            fi("TU", "TU_1_2", 1);
            fi("TU", "TU_2_1", 2); // Tumor Tier-1
            fi("TU", "TU_2_2", 3);
            
            /*
             * Eg: "AD_1_1" -> first value in the first sample
             *     "AD_2_1" -> firsr value in the second sample
             *
             * Note that we assume the sample is diploid.
             */
            
            fi("AD", "AD_1_1", 0);
            fi("AD", "AD_1_2", 1);
            fi("AD", "AD_2_1", 2);
            fi("AD", "AD_2_2", 3);
            
            x.hdr  = (void *) hdr;
            x.line = (void *) line;
            
            if (x.ifs.count("GT"))
            {
                if      (x.ifs["GT"] == "MSI") { x.gt = Genotype::MSI; }
                else if (x.ifs["GT"] == "HOM") { x.gt = Genotype::Homozygous;  }
                else if (x.ifs["GT"] == "HET") { x.gt = Genotype::Heterzygous; }
                else if (x.ifs["GT"] == "SOM") { x.gt = Genotype::Somatic;     }
            }
            
            /*
             * Depth metrics
             */
            
            int32_t *a1 = NULL, a2 = 0;
            if (bcf_get_format_int32(hdr, line, "AD", &a1, &a2) >= 2)
            {
                for (auto i = 0; i < a2 && i < 4; i++) { x.AD[i] = a1[i]; }
                free(a1);
            }
            
            /*
             * Quality metrics
             */
            
            x.qual[0] = line->qual;
            
            if (isStrelka(x) && x.iff.count("SomaticEVS"))
            {
                /*
                 * Strelka has the QSI, QSS and SomaticEVS metrics. SomaticEVS is choosen
                 * because it's defined for both SNPs and indels. Set both normal and tumor
                 * quality score equal.
                 */
                
                x.qual[0] = x.qual[1] = x.iff.at("SomaticEVS");
            }
            else if (x.iff.count("TLOD"))
            {
                x.qual[0] = x.qual[1] = x.iff.at("TLOD");
            }

            f(x);
        };
        
        /*
        if (line->n_allele <= 2)
        {
            for (auto i = 1; i < line->n_allele; i++)
            {
                run(i);
            }
        }
         */

        run(1);
    }
    
    free(ps);
    free(pi);
    free(pf);
    hts_close(fp);
}
