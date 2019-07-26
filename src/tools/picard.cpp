#include <sstream>
#include "tools/picard.hpp"

using namespace Anaquin;

std::string Picard::report() const
{
    std::stringstream ss;
    
    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\n";
    ss << ((boost::format(f) % "CHROM" % "POS" % "REF" % "MATCHES" % "MISMATCHES" % "DELETIONS" % "INSERTIONS" % "A" % "C" % "G" % "T").str());
    
    // For each reference genome ...
    for (auto &i : seqs)
    {
        const auto &cID = i.first;
        
        // For each position in the reference genome ...
        for (auto j = 0; j < i.second.size(); j++)
        {
            Counts A = 0; // Alternative allele
            Counts C = 0; // Alternative allele
            Counts G = 0; // Alternative allele
            Counts T = 0; // Alternative allele
            
            switch (i.second[j])
            {
                case 'a':
                case 'A':
                {
                    C = snps[cID][j][SNPBin::AC];
                    G = snps[cID][j][SNPBin::AG];
                    T = snps[cID][j][SNPBin::AT];
                    break;
                }
                    
                case 'c':
                case 'C':
                {
                    A = snps[cID][j][SNPBin::CA];
                    T = snps[cID][j][SNPBin::CT];
                    G = snps[cID][j][SNPBin::CG];
                    break;
                }
                    
                case 'g':
                case 'G':
                {
                    A = snps[cID][j][SNPBin::GC];
                    C = snps[cID][j][SNPBin::GC];
                    T = snps[cID][j][SNPBin::GT];
                    break;
                }
                    
                case 't':
                case 'T':
                {
                    A = snps[cID][j][SNPBin::TC];
                    C = snps[cID][j][SNPBin::TC];
                    G = snps[cID][j][SNPBin::TG];
                    break;
                }
                    
                default : { continue; }
            }
            
            const auto YY = snps[cID][j][SNPBin::Match];
            const auto NN = A + C + G + T;
            
            const auto D = sum(dls[cID][j]);
            const auto I = sum(ins[cID][j]);
            
            ss << ((boost::format(f) % cID % j % i.second[j] % YY  % NN % D % I % A % C % G % T).str());
        }
    }
    
    return ss.str();
}
