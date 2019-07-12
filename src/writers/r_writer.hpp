#ifndef R_WRITER_HPP
#define R_WRITER_HPP

#include <set>
#include <math.h>
#include <numeric>
#include "stats/linear.hpp"
#include <boost/format.hpp>

// Defined in main.cpp
extern std::string __full_command__;

namespace Anaquin
{
    struct RWriter
    {
        static Scripts createGene(const FileName    &,
                                  const Path        &,
                                  const std::string &,
                                  const std::string &,
                                  const std::string &,
                                  const std::string &,
                                  const std::string &,
                                  bool showLOQ = false,
                                  bool shouldLog = true);

        static Scripts createLinear(const FileName    &,
                                    const Path        &,
                                    const std::string &,
                                    const std::string &,
                                    const std::string &,
                                    const std::string &,
                                    const std::string &,
                                    bool showLOQ = false,
                                    bool shouldLog = true,
                                    const std::string &script = "");

        static Scripts createScript(const FileName &, const Scripts &);
        static Scripts createScript(const FileName &, const Scripts &, const std::string &);
    };
}

#endif
