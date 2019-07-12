#ifndef BAM_WRITER_HPP
#define BAM_WRITER_HPP

#include "data/data.hpp"
#include "parsers/parser_bam.hpp"

namespace Anaquin
{
    class BAMWriter
    {
        public:        
            BAMWriter();

            void close();
            void open(const FileName &);
        
            void write (const ParserBAM::Data &);
            void writeH(const ParserBAM::Data &);

        private:
            struct Impl;
            Impl *_impl;
    };
}

#endif
