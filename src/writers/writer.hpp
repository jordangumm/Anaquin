#ifndef WRITER_HPP
#define WRITER_HPP

#include "data/data.hpp"

namespace Anaquin
{
    template <typename T = std::string> struct Writer
    {
        virtual void close() = 0;
        virtual void open(const FileName &) = 0;
        virtual void write(const T &) = 0;
    };
    
    struct MockWriter : public Writer<>
    {
        inline void close() override {}
        inline void open(const FileName &) override {}
        inline void write(const std::string &) override {}
    };
}

#endif
