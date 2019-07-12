#include "tools/tools.hpp"
#include "data/resources.hpp"
#include "writers/r_writer.hpp"

using namespace Anaquin;

// Defined in main.cpp
extern Path __output__;

Scripts RWriter::createGene(const FileName    &file,
                            const Path        &path,
                            const std::string &title,
                            const std::string &xl,
                            const std::string &yl,
                            const std::string &exp,
                            const std::string &obs,
                            bool  showLOQ,
                            bool  shouldLog)
{
    return (boost::format(PlotGene())
                                % date()
                                % __full_command__
                                % path
                                % file
                                % title
                                % xl
                                % yl
                                % exp
                                % obs
                                % (showLOQ ? "T" : "F")
                                % (shouldLog ? "T" : "F")).str();
}

Scripts RWriter::createLinear(const FileName    &file,
                              const Path        &path,
                              const std::string &title,
                              const std::string &xl,
                              const std::string &yl,
                              const std::string &exp,
                              const std::string &obs,
                              bool  showLOQ,
                              bool  shouldLog,
                              const std::string &script)
{
    return (boost::format(script.empty() ? PlotLinear() : script)
                                  % date()
                                  % __full_command__
                                  % path
                                  % file
                                  % title
                                  % xl
                                  % yl
                                  % exp
                                  % obs
                                  % (showLOQ ? "T" : "F")
                                  % (shouldLog ? "T" : "F")).str();
}

Scripts RWriter::createScript(const FileName &file, const Scripts &script)
{
    return (boost::format(script) % date()
                                  % __full_command__
                                  % __output__
                                  % file).str();
}

Scripts RWriter::createScript(const FileName &file, const Scripts &script, const std::string &x)
{
    return (boost::format(script) % date()
                                  % __full_command__
                                  % __output__
                                  % file
                                  % x).str();
}
