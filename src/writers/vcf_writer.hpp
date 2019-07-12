#ifndef VCF_WRITER_HPP
#define VCF_WRITER_HPP

#include "data/data.hpp"
#include "writers/writer.hpp"

namespace Anaquin
{
    class VCFWriter : public Writer<>
    {
        public:
        
            ~VCFWriter() { close(); }

            void open(const FileName &) override;

            void write(void *hdr, void *line);
            void write(const std::string &) override {}

        private:

            void close() override;

            bool _head;
            void *_fp;
    };
}

#endif
