#include <vector>
#include "data/fastq.hpp"
#include "tools/tools.hpp"

using namespace Anaquin;

struct FASTQ::Impl
{
    std::vector<std::string> heads;
};

FASTQ::FASTQ()
{
    _impl = std::shared_ptr<FASTQ::Impl>(new FASTQ::Impl());
}

unsigned FASTQ::heads() const { return _impl->heads.size(); }

std::string FASTQ::inst(Format f) const
{
    switch (f)
    {
        case Illumina_V2: { return tokN(_impl->heads.front(), ":", 0); }
        default : { return "-"; }
    }
}

std::string FASTQ::run(Format f) const
{
    switch (f)
    {
        case Illumina_V2: { return tokN(_impl->heads.front(), ":", 1); }
        default : { return "-"; }
    }
}

std::string FASTQ::flow(Format f) const
{
    switch (f)
    {
        case Illumina_V2: { return tokN(_impl->heads.front(), ":", 2); }
        default : { return "-"; }
    }
}

std::string FASTQ::lane(Format f) const
{
    switch (f)
    {
        case Illumina_V2: { return tokN(_impl->heads.front(), ":", 3); }
        default : { return "-"; }
    }
}

FASTQ::Format FASTQ::format()
{
    if (_impl->heads.empty())
    {
        return Format::None;
    }

    auto tmp1 = std::vector<std::string>();
    split(_impl->heads.front(), " ", tmp1);

    auto tmp2 = std::vector<std::string>();
    split(tmp1.front(), ":", tmp2);
    
    if (tmp2.size() == 7)
    {
        return Format::Illumina_V2;
    }

    return Format::None;
}

void FASTQ::addHead(const std::string &x)
{
    _impl->heads.push_back(x);
}
