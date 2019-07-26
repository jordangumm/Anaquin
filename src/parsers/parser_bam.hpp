#ifndef PARSER_BAM_HPP
#define PARSER_BAM_HPP

#include <memory>
#include "data/alignment.hpp"
#include "stats/analyzer.hpp"

namespace Anaquin
{
    struct ParserBAM
    {
        struct Info
        {
            Progress p = 0;

            // Whether this is a multi-alignment
            bool multi;
            
            // Whether there is insertion
            bool ins;
            
            // Whether there is deletion
            bool del;
            
            // Whether there is skipped region
            bool skip;
            
            // Whether there is clipping
            bool clip;
            
            // Size of the chromosome of the alignment
            Base length;
            
            void *b;
            void *h;
        };

        class Data : public Alignment
        {
            friend struct ParserBAM;
            
            public:
            
                bool nextCigar(Locus &l, bool &spliced);

                /*
                 * Optional fields
                 */
            
                void lSeq();     // Lazy loadding of sequence
                void lQual();    // Lazy loading of quality
                void lName();    // Lazy loading of reads
                void lMateID();  // Lazy loading of mate ID
                void lMatePos(); // Lazy loading of mate position
                void lCigar();

                inline void *b() const { return _b; } // bam1_t
                inline void *h() const { return _h; } // bam_hdr_t

                void *copyH() const;
                void *copyB() const;

            private:
            
                mutable unsigned _i, _n;

                void *_b;
                void *_h;
        };
        
        typedef std::function<void (Data &, const Info &)> Functor;
        
        static std::map<ChrID, Base> header(const FileName &);
        
        /*
         * In order to improve the efficiency, not everything is computed. Set the last
         * argument to true will force it to happen.
         */

        static void parse(const FileName &, Functor);
        
        struct ParseResult
        {
            // Sequin read?
            bool isSeq = false;
            
            /*
             * Should we reverse complement this read sequence and
             * also reverse the quality scores?
             *
             *     std::reverse(x.qual.begin(), x.qual.end());
             *     x.seq = revcomp(x.seq);
             */
            
            bool rc = false;
        };
        
        // Functor for parsing alignments selectively
        typedef std::function<ParseResult (Data &, const Info &)> F2;

        static void parse_(const FileName &, AlignedReads &r, F2);
    };
}

#endif
