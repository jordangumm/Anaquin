#ifndef FASTQ_HPP
#define FASTQ_HPP

#include <string>
#include <memory>

namespace Anaquin
{
    struct FASTQ
    {
        enum Format
        {
            None,
            Illumina_V2
        };
        
        struct Impl;
        FASTQ();
        
        // Add a FASTQ header
        void addHead(const std::string &);
        
        // Number of headers already added
        unsigned heads() const;

        Format format();
        
        std::string run(Format)  const;
        std::string inst(Format) const;
        std::string flow(Format) const;
        std::string lane(Format) const;

        std::shared_ptr<Impl> _impl;
    };
}

#endif
