#ifndef SAM_WRITER_HPP
#define SAM_WRITER_HPP

#include <htslib/sam.h>
#include "tools/samtools.hpp"
#include "writers/writer.hpp"
#include "parsers/parser_bam.hpp"

namespace Anaquin
{
    class SAMWriter : public Writer<ParserBAM::Data>
    {
        public:

            SAMWriter(bool cout, bool head) : _cout(cout), _head(head) {}
        
            inline void close() override
            {
                sam_close(_fp);
            }

            inline void open(const FileName &file) override
            {
                if (_cout) { _fp = sam_open("/dev/null",  "w"); }
                else       { _fp = sam_open(file.c_str(), "w"); }
            }

            inline void write(const ParserBAM::Data &x) override
            {
                const auto *b = reinterpret_cast<bam1_t *>(x.b());
                const auto *h = reinterpret_cast<bam_hdr_t *>(x.h());

                if (_head)
                {
                    if (_cout) { std::cout << std::string(h->text); }
                    else
                    {
                        if (sam_hdr_write(_fp, h) == -1)
                        {
                            throw std::runtime_error("sam_hdr_write() failed");
                        }
                    }
                }
                
                if (sam_write1(_fp, h, b) == -1)
                {
                    throw std::runtime_error("Failed to SAM record");
                }
                
#ifndef DEBUG
                if (_cout) { std::cout << std::string(_fp->line.s); }
#endif
                _head = false;
            }

        private:

            bool _cout = true;
        
            // Write header?
            bool _head = false;

            // File pointer
            samFile *_fp;
    };
}

#endif
