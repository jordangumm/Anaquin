#include "tools/tools.hpp"
#include "tools/report.hpp"
#include "data/resources.hpp"

using namespace Anaquin;

void Report::gReport(const Options &o)
{
    A_ASSERT(o.v);
    o.generate("germline_report.html");
    runScript(PReport(), ((boost::format("germline %1% %2% %3% %4% %5%") % o.work
                                                                         % o.work
                                                                         % "germline"
                                                                         % script2File(GermHTML())
                                                                         % script2File(CSS()))).str());
}

void Report::sReport(const Options &o)
{
    A_ASSERT(o.v);
    o.generate("somatic_report.html");
    runScript(PReport(), ((boost::format("somatic %1% %2% %3% %4% %5%") % o.work
                                                                        % o.work
                                                                        % "somatic"
                                                                        % script2File(SomaHTML())
                                                                        % script2File(CSS()))).str());
}

void Report::calibrate(const Options &o)
{
    A_ASSERT(o.k);
    o.generate("calibrate_report.html");
    runScript(PReport(), ((boost::format("calibrate %1% %2% %3% %4% %5%") % o.work
                                                                          % o.work
                                                                          % "calibrate"
                                                                          % script2File(GSplitHTML())
                                                                          % script2File(CSS()))).str());
}

void Report::genome(const Options &o)
{
    A_ASSERT(o.k);
    o.generate("split_report.html");
    runScript(PReport(), ((boost::format("genome %1% %2% %3% %4% %5%") % o.work
                                                                       % o.work
                                                                       % "split"
                                                                       % script2File(GSplitHTML())
                                                                       % script2File(CSS()))).str());
}

void Report::rna(const Options &o)
{
    A_ASSERT(o.k);
    o.generate("rna_report.html");
    runScript(PReport(), ((boost::format("rna %1% %2% %3% %4% %5%") % o.work
                                                                    % o.work
                                                                    % "rna"
                                                                    % script2File(RSplitHTML())
                                                                    % script2File(CSS()))).str());
}

void Report::meta(const Options &o)
{
    A_ASSERT(o.k);
    o.generate("meta_report.html");
    runScript(PReport(), ((boost::format("meta %1% %2% %3% %4% %5%") % o.work
                                                                     % o.work
                                                                     % "meta"
                                                                     % script2File(MSplitHTML())
                                                                     % script2File(CSS()))).str());
}
