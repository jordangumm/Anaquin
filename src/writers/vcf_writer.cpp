#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "writers/vcf_writer.hpp"

using namespace Anaquin;

void VCFWriter::open(const FileName &file)
{
    _fp = bcf_open(file.c_str(), "w");
}

void VCFWriter::write(void *hdr, void *line)
{
    if (!_head)
    {
        _head = true;
        bcf_hdr_write((htsFile *) _fp, (bcf_hdr_t *) hdr);
    }

    bcf_write1((htsFile *) _fp, (bcf_hdr_t *) hdr, (bcf1_t *) line);
}

void VCFWriter::close()
{
    if (_fp) { bcf_close((htsFile *) _fp); }
}
