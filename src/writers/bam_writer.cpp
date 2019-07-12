#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include "writers/bam_writer.hpp"

using namespace Anaquin;

struct BAMWriter::Impl
{
    BGZF *f;
    bool head; // Should we write header?
};

BAMWriter::BAMWriter()
{
    _impl = new Impl();
}

void BAMWriter::close()
{
    bgzf_close(_impl->f);
}

void BAMWriter::open(const FileName &file)
{
    _impl->f = bgzf_open(file.data(), "w");
}

void BAMWriter::writeH(const ParserBAM::Data &x)
{
    if (bam_hdr_write(_impl->f, reinterpret_cast<bam_hdr_t *>(x.h())) == -1)
    {
        throw std::runtime_error("bam_hdr_write() failed");
    }
    
    _impl->head = true;
}

void BAMWriter::write(const ParserBAM::Data &x)
{
    if (!_impl->head && bam_hdr_write(_impl->f, reinterpret_cast<bam_hdr_t *>(x.h())) == -1)
    {
        throw std::runtime_error("bam_hdr_write() failed");
    }
    
    _impl->head = true;
    
    if (bam_write1(_impl->f, (bam1_t *) x.b()) == -1)
    {
        throw std::runtime_error("bam_write1() failed");
    }
}
