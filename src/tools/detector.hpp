#ifndef DETECTOR_HPP
#define DETECTOR_HPP

#include "data/data.hpp"

namespace Anaquin
{
    struct Detector
    {
        static HumanAssembly fromBAM(const FileName &);
        static HumanAssembly fromVCF(const FileName &);
    };
}

#endif
