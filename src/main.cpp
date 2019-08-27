#include <ctime>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <strings.h>
#include <execinfo.h>
#include <sys/stat.h>

#include "RNA/RNA.hpp"
#include "RNA/r_split.hpp"

#include "Meta/Meta.hpp"
#include "Meta/m_split.hpp"

#include "Genomics/g_germ.hpp"
#include "Genomics/g_split.hpp"
#include "Genomics/g_somatic.hpp"
#include "Genomics/g_calibrate.hpp"

#include "parsers/parser_csv.hpp"
#include "parsers/parser_vcf.hpp"

#include "tools/tools.hpp"
#include "tools/report.hpp"
#include "data/resources.hpp"
#include "tools/detector.hpp"
#include "tools/bedtools.hpp"
#include "tools/checkPackage.h"
#include "tools/attributes.hpp"
#include "writers/file_writer.hpp"
#include "writers/terminal_writer.hpp"

#ifdef UNIT_TEST
#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#endif

#define DEFAULT_EDGE 550

using namespace Anaquin;

typedef int Option;

typedef std::string Value;

static HumanAssembly __ha__;

static std::string version() { return "3.9.0"; }

/*
 * Options specified in the command line
 */

#define OPT_TEST    320
#define OPT_TOOL    321
#define OPT_PATH    325

#define OPT_R_BED    801
#define OPT_METHOD   802
#define OPT_R_VCF    803
#define OPT_COMBINE  804
#define OPT_TRIM     805
#define OPT_MIXTURE  806
#define OPT_R_HUMAN  807
#define OPT_R_DECOY  808
#define OPT_R_REGS   809
#define OPT_RESOURCE 810
#define OPT_BUILD    811
#define OPT_SYNC     812
#define OPT_ATTTSV   813
#define OPT_FASTA    814
#define OPT_KMER     815
#define OPT_RULE     816
#define OPT_THREAD   817
#define OPT_EDGE     818
#define OPT_SKIP     819
#define OPT_MERGE    820
#define OPT_WINDOW   821
#define OPT_1        822
#define OPT_2        823
#define OPT_REPORT   824
#define OPT_SAMPLE   825
#define OPT_SEQUIN   826
#define OPT_S_FRAC   827
#define OPT_S_ABSO   828
#define OPT_L_FRAC   829
#define OPT_L_ABSO   830
#define OPT_FASTQ    831
#define OPT_NO_FLIP  832
#define OPT_MMIX     833

// Shared with other modules
std::string __full_command__;

static Path __working__;

// Shared with other modules
Path __output__;

std::shared_ptr<FileWriter> __loggerW__;
std::shared_ptr<TerminalWriter> __outputW__;

static std::map<Value, Tool> _tools =
{
    { "rna",       Tool::RNA       },
    { "meta",      Tool::Meta      },
    { "split",     Tool::Split     },
    { "germline",  Tool::Germline  },
    { "somatic",   Tool::Somatic   },
    { "calibrate", Tool::Calibrate },
};

struct Parsing
{
    Path path = "output";

    // Specific options
    std::map<Option, std::string> opts;
    
    std::map<Option, double> od;
    
    // How Anaquin is invoked
    Command  cmd;

    // Mixture A(1) or mixture B(2)?
    Mixture mix = Mix_1;
    
    Tool tool;
};

// Wrap the variables so that it'll be easier to reset them
static Parsing _p;

static Scripts fixManual(const Scripts &str)
{
    auto x = str;
    
    boost::replace_all(x, "<b>", "\e[1m");
    boost::replace_all(x, "</b>", "\e[0m");
    boost::replace_all(x, "<i>", "\e[3m");
    boost::replace_all(x, "</i>", "\e[0m");
    
    return x;
}

struct InvalidUsageException : public std::exception {};

struct InvalidDependency : public std::runtime_error
{
    InvalidDependency(const std::string &x) : std::runtime_error(x) {}
};

struct InvalidOptionException : public std::exception
{
    InvalidOptionException(const std::string &o) : o(o) {}

    const std::string o;
};

struct InvalidValueException : public std::exception
{
    InvalidValueException(const std::string &o, const std::string &v) : o(o), v(v) {}

    const std::string o, v;
};

struct InvalidToolError : public InvalidValueException
{
    InvalidToolError(const std::string &v) : InvalidValueException("-t", v) {}
};

struct MissingOptionError : public std::exception
{
    MissingOptionError(const std::string &o) : o(o) {}

    // Option that is missing
    const std::string o;
    
    // Possible values for the missing option
    const std::string r;
};

struct MissingBundleError : public std::runtime_error
{
    MissingBundleError(const Path &path) : std::runtime_error(path), path(path) {}
    const Path path;
};

struct InvalidBundleError : public std::runtime_error
{
    InvalidBundleError(const FileName &file) : std::runtime_error(file), file(file) {}
    const FileName file;
};

static const char *short_opts = ":";

static const struct option long_opts[] =
{
    { "1", required_argument, 0, OPT_1 },
    { "2", required_argument, 0, OPT_2 },

    { "t",       required_argument, 0, OPT_THREAD },
    { "threads", required_argument, 0, OPT_THREAD },

    { "sample",   required_argument, 0, OPT_SAMPLE  },
    { "sequin",   required_argument, 0, OPT_SEQUIN  },
    { "combined", required_argument, 0, OPT_COMBINE },
    
    { "fq",      no_argument, 0, OPT_FASTQ   },
    { "merge",   no_argument, 0, OPT_MERGE   },
    { "report",  no_argument, 0, OPT_REPORT  },
    { "no_flip", no_argument, 0, OPT_NO_FLIP },
    
    { "noflip",  no_argument, 0, OPT_NO_FLIP  },
    { "human_regions",    required_argument, 0, OPT_R_HUMAN },
    { "decoy_regions",    required_argument, 0, OPT_R_DECOY },
    { "restrict_regions", required_argument, 0, OPT_R_REGS  },

    { "manual_mix",     required_argument, 0, OPT_MMIX  },
    { "manual_bed",     required_argument, 0, OPT_R_BED },
    { "decoy_variants", required_argument, 0, OPT_R_VCF },
    { "manual_fasta",   required_argument, 0, OPT_FASTA },

    { "kmer",      required_argument, 0, OPT_KMER   },
    { "threshold", required_argument, 0, OPT_RULE   },
    { "skip",      required_argument, 0, OPT_SKIP   },

    { "r",            required_argument, 0, OPT_RESOURCE },
    { "resource_dir", required_argument, 0, OPT_RESOURCE },

    { "build", required_argument, 0, OPT_BUILD }, // Alternative genome assembly

    { "mix",    required_argument, 0, OPT_MIXTURE },
    { "trim",   required_argument, 0, OPT_TRIM    },
    { "method", required_argument, 0, OPT_METHOD  },
    
    { "sequin_fraction", required_argument, 0, OPT_S_FRAC },
    { "sequin_absolute", required_argument, 0, OPT_S_ABSO },
    { "ladder_fraction", required_argument, 0, OPT_L_FRAC },
    { "ladder_absolute", required_argument, 0, OPT_L_ABSO },
    
    { "edge",   required_argument, 0, OPT_EDGE   },
    { "window", required_argument, 0, OPT_WINDOW },

    { "o",      required_argument, 0, OPT_PATH },
    { "output", required_argument, 0, OPT_PATH },

    {0, 0, 0, 0 }
};

static HumanAssembly parseBuild(const std::string &x)
{
    if      (x == "hg38") { return HumanAssembly::hg38; }
    else if (x == "hg19") { return HumanAssembly::hg19; }
    else if (x == "gr37") { return HumanAssembly::gr37; }
    else if (x == "gr38") { return HumanAssembly::gr38; }
    else
    {
        throw InvalidValueException("--build", "Invalid assembly build. Must be \"hg38\", \"hg19\", \"gr37\" or \"gr38\"");
    }
}

static std::string optToStr(int opt)
{
    for (const auto &o : long_opts)
    {
        if (o.val == opt)
        {
            return o.name;
        }
    }
    
    throw std::runtime_error("Invalid option: " + std::to_string(opt));
}

/*
static void checkReport()
{
    if (!PackageChecker::checkR())
    {
        throw InvalidDependency("R installation is required for report generation. Please install R on https://www.r-project.org.");
    }
    else if (!PackageChecker::checkPython())
    {
        throw InvalidDependency("Python installation is required for report generation. Please install Python on https://www.python.org.");
    }
    else if (!PackageChecker::checkRAnaquin(version()))
    {
        throw InvalidDependency("Bioconductor Anaquin R-package version " + version() + " is required for report generation.");
    }
}
*/

static void printUsage()
{
    extern Scripts Manual();
    std::cout << std::endl << fixManual(Manual()) << std::endl << std::endl;
}

static Scripts manual(Tool tool)
{
    switch (tool)
    {
        case Tool::RNA:       { return rna();       }
        case Tool::Meta:      { return meta();      }
        case Tool::Split:     { return split();     }
        case Tool::Somatic:   { return somatic();   }
        case Tool::Germline:  { return germline();  }
        case Tool::Calibrate: { return calibrate(); }
        default:              { return "";          }
    }
}

inline std::string option(const Option &key, const Scripts &x = "")
{
    return _p.opts.count(key) ? _p.opts[key] : x;
}

template <typename F> std::shared_ptr<Ladder> readL(F f, Option key, UserReference &, const Scripts &x = "")
{
    return _p.opts.count(key) ? std::shared_ptr<Ladder>(new Ladder(f(Reader(_p.opts[key])))) :
                                std::shared_ptr<Ladder>(new Ladder(f(Reader(x))));
}

static std::shared_ptr<VCFLadder> readV(Option opt, UserReference &, std::shared_ptr<BedData> rb, const FileName &def)
{
    auto rr = (_p.opts.count(opt)) ? Reader(_p.opts[opt]) : Reader(def);
    return std::shared_ptr<VCFLadder>(new VCFLadder(Standard::addVCF(rr, rb)));
}

static std::shared_ptr<VCFLadder> readGV(Option opt, UserReference &, std::shared_ptr<BedData> rb, const FileName &def)
{
    auto rr = (_p.opts.count(opt)) ? Reader(_p.opts[opt]) : Reader(def);
    return std::shared_ptr<VCFLadder>(new VCFLadder(Standard::addGVCF(rr, rb)));
}

static std::shared_ptr<VCFLadder> readSV(Option opt, UserReference &, std::shared_ptr<BedData> rb, const FileName &def)
{
    auto rr = (_p.opts.count(opt)) ? Reader(_p.opts[opt]) : Reader(def);
    return std::shared_ptr<VCFLadder>(new VCFLadder(Standard::addSVCF(rr, rb)));
}

static std::shared_ptr<BedData> readR(const FileName &file, const RegionOptions &o = RegionOptions())
{
    return std::shared_ptr<BedData>(new BedData(Standard::readBED(Reader(file), o)));
}

static WriterOptions *__o__;

void info(const std::string &s) { if (__o__) { __o__->info(s); } }
void wait(const std::string &s) { if (__o__) { __o__->wait(s); } }

template <typename Analyzer, typename O, typename F> void start(const std::string &name, F f, O &o)
{
    o.info("-----------------------------------------");
    o.info("------------- Sequin Analysis -----------");
    o.info("-----------------------------------------");

    const auto path = _p.path;

    // Required for logging
    __o__ = &o;

    o.output = __outputW__;
    o.logger = __loggerW__;
    o.writer = std::shared_ptr<FileWriter>(new FileWriter(path));
    
    if (system(("mkdir -p " + path).c_str()))
    {
        throw std::runtime_error("Failed to create output directory");
    }
    
    o.name = name;
    o.cmd = __full_command__ = _p.cmd;
    o.work = path;
    
    assert(!o.cmd.empty());
    assert(!o.name.empty());
    assert(!o.work.empty());

    o.info("Version: " + version());
    o.info(_p.cmd);
    o.info(date());
    o.info("Path: " + path);
    o.info("Resources: " + _p.opts[OPT_RESOURCE]);

    if (o.report)
    {
        createD(o.work + "/report_files");
    }
    
    using namespace std::chrono;
    
    auto begin = high_resolution_clock::now();

    f(o);
    
    auto end = high_resolution_clock::now();

    // Remove all working directories
    clearAllTmp();
    
    const auto elapsed = (boost::format("Completed. %1% seconds.") % duration_cast<seconds>(end - begin).count()).str();
    o.info(elapsed);

    o.logger->close();
}

template <typename Analyzer> void analyze_1(const std::string &p, Option x, typename Analyzer::Options o = typename Analyzer::Options())
{
    return start<Analyzer>(p, [&](const typename Analyzer::Options &o)
    {
        Analyzer::report(_p.opts.at(x), o);
    }, o);
}

template <typename Analyzer> void analyze_2(const std::string &p, Option x1, Option x2, typename Analyzer::Options &o = typename Analyzer::Options())
{
    return start<Analyzer>(p, [&](const typename Analyzer::Options &o)
    {
        Analyzer::report(_p.opts.count(x1) ? _p.opts[x1] : "", _p.opts.count(x2) ? _p.opts[x2] : "", o);
    }, o);
}

static std::shared_ptr<Translation> readTrans(const Reader &r)
{
    Translation l = Translation();
    
    ParserCSV::parse(r, [&](const ParserCSV::Data &x, Progress)
    {
        if (x.size() < 2)
        {
            throw std::runtime_error("Invalid format. Two or more columns expected");
        }
        
        l[x[0]] = x[1];
    }, "\t");

    return std::shared_ptr<Translation>(new Translation(l));
}

template <typename T> std::shared_ptr<Ladder> readTSV(const Reader &r, T &, int con = 1)
{
    Ladder l = Ladder();
    
    ParserCSV::parse(r, [&](const ParserCSV::Data &x, Progress)
    {
        if (x.size() < 2)
        {
            throw std::runtime_error("Invalid format. Two or more columns expected");
        }
        else if (x[con] == MISSING || !isNumber(x[con]))
        {
            return;
        }
        
        l.add(x[0], Mix_1, stof(x[con]));
    }, "\t");

    return std::shared_ptr<Ladder>(new Ladder(l));
}

void parse(int argc, char ** argv)
{
    auto tmp = new char*[argc+1];
    
    for (auto i = 0; i < argc; i++)
    {
        tmp[i] = (char *) malloc((strlen(argv[i]) + 1) * sizeof(char));
        strcpy(tmp[i], argv[i]);
    }
    
    _p = Parsing();

    if ((argc <= 1) || (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")))
    {
        printUsage();
        return;
    }

    // Required for unit-testing
    optind = 1;

    /*
     * Reconstruct the overall command
     */
    
    for (auto i = 0; i < argc; i++)
    {
        _p.cmd += std::string(argv[i]) + " ";
    }

    A_ASSERT(!_p.cmd.empty());

    int next, index;

    auto checkPath = [&](const Path &path)
    {
        if (path[0] == '/')
        {
            return path;
        }
        else
        {
            return __working__ + "/" + path;
        }
    };
    
    auto checkBundle = [&](const Path &path)
    {
        if (!exists(path))
        {
            throw MissingBundleError(path);
        }
        
        const auto files = listFiles(path);
        
        auto check = [&](const FileName &x)
        {
            if (std::find_if(files.begin(), files.end(), [&](const FileName &i) {
                return isBegin(i, x);
            }) == files.end()) { throw InvalidBundleError(x); }
        };
        
        check("genome");
        check("synthetic");
        check("metagenome");
        check("transcriptome");

        return path;
    };

    auto checkFile = [&](const FileName &file)
    {
        if (!std::ifstream(file).good())
        {
            throw InvalidFileError(file);
        }
        
        return file;
    };

    /*
     * Pre-process arguments. This way, we can examine the options in whatever order we'd like to
     */

    std::vector<Value>  vals;
    std::vector<Option> opts;

    if (strcmp(argv[1], "-v") == 0)
    {
        std::cout << version() << std::endl;
        return;
    }
    else if (strcmp(argv[1], "-t") == 0)
    {
#ifdef UNIT_TEST
        Catch::Session().run(1, argv);
#else
        A_THROW("UNIT_TEST is undefined");
#endif
        return;
    }
    else if (!_tools.count(argv[1]))
    {
        throw InvalidToolError(argv[1]);
    }

    _p.tool = _tools[argv[1]];
    const auto isHelp = argc >= 3 && (!strcmp(argv[2], "-h") || !strcmp(argv[2], "--help"));

    if (isHelp)
    {
        if (argc != 3)
        {
            throw std::runtime_error("Too many arguments for help usage. Usage: anaquin <tool> -h or anaquin <tool> --help");
        }

        std::cout << std::endl << fixManual(manual(_p.tool)) << std::endl << std::endl;
        return;
    }
    
    unsigned n = 2;

    while ((next = getopt_long_only(argc, argv, short_opts, long_opts, &index)) != -1)
    {
        if (next < OPT_TOOL)
        {
            throw InvalidOptionException(argv[n]);
        }
        
        opts.push_back(next);
        
        // Whether this option has an value
        const auto hasValue = optarg;
        
        n += hasValue ? 2 : 1;
        
        vals.push_back(hasValue ? std::string(optarg) : "");
    }

    for (auto i = 0u; i < opts.size(); i++)
    {
        auto opt = opts[i];
        auto val = vals[i];

        switch (opt)
        {
            case OPT_RULE:
            {
                try
                {
                    stof(val);
                    _p.opts[opt] = val;
                }
                catch (...)
                {
                    throw std::runtime_error(val + " is not an. Please check and try again.");
                }
                
                break;
            }

            case OPT_S_FRAC:
            case OPT_S_ABSO:
            case OPT_L_FRAC:
            case OPT_L_ABSO:
            {
                try
                {
                    _p.od[opt] = stod(val);
                }
                catch (...)
                {
                    throw std::runtime_error(val + " is not a floating number");
                }

                break;
            }
                
            case OPT_EDGE:
            case OPT_KMER:
            case OPT_SKIP:
			case OPT_WINDOW:
            {
                try
                {
                    stoi(val);
                    _p.opts[opt] = val;
                }
                catch (...)
                {
                    throw std::runtime_error(val + " is not an integer. Please check and try again.");
                }

                break;
            }

            case OPT_METHOD:
            {
                switch (_p.tool)
                {
                    case Tool::Somatic:
                    case Tool::Germline:
                    case Tool::Calibrate: { _p.opts[opt] = val; break; }
                    default : { throw InvalidOptionException("Invalid usage for --method"); }
                }
                
                break;
            }

            case OPT_TRIM:
            case OPT_FASTQ:
            case OPT_FASTA:
            case OPT_BUILD:
            case OPT_MERGE:
            case OPT_THREAD:
            case OPT_REPORT:
            case OPT_NO_FLIP:
            {
                if (opt == OPT_REPORT)
                {
                    //checkReport();
                }
                
                _p.opts[opt] = val;
                break;
            }

            case OPT_MIXTURE:
            {
                if      (val == "A") { _p.mix = Mixture::Mix_1; }
                else if (val == "B") { _p.mix = Mixture::Mix_2; }
                else if (val == "C") { _p.mix = Mixture::Mix_3; }
                else                 { throw InvalidValueException("-mix", val); }
                break;
            }

            case OPT_1:
            case OPT_2:
            case OPT_SAMPLE:
            case OPT_SEQUIN:
            case OPT_COMBINE: { checkFile(_p.opts[opt] = val); break; }

            case OPT_RESOURCE: { checkBundle(_p.opts[opt] = val); break; }
                
            case OPT_MMIX:
            case OPT_R_VCF:
            case OPT_R_BED:
            case OPT_R_REGS:
            case OPT_R_HUMAN:
            case OPT_R_DECOY:
            {
                checkFile(_p.opts[opt] = val);
                break;
            }

            case OPT_PATH: { _p.path = val; break; }

            default: { throw InvalidUsageException(); }
        }
    }
    
    if (!_p.opts.count(OPT_RESOURCE))
    {
        _p.opts[OPT_RESOURCE] = checkBundle(execPath() + "/resources");
    }

    __output__ = _p.path = checkPath(_p.path);
    
    if (opts.empty())
    {
        std::cout << std::endl << fixManual(manual(_p.tool)) << std::endl << std::endl;
        return;
    }
    
    UserReference r;
    
    auto &s = Standard::instance();
    
    __loggerW__ = std::shared_ptr<FileWriter>(new FileWriter(_p.path));
    __outputW__ = std::shared_ptr<TerminalWriter>(new TerminalWriter());
    __loggerW__->open("anaquin.log");
    
    #define GSeqFA_()      (GSeqFA(_p.opts.at(OPT_RESOURCE)).path)
    #define GFeatBED_(x)   (GFeatBED(_p.opts.at(OPT_RESOURCE), x).path)
    #define GSynTSV_(x)    (GSynTSV(_p.opts.at(OPT_RESOURCE)).path)
    #define GAttrBED_(x)   (GAttrBED(_p.opts.at(OPT_RESOURCE)).path)
    #define GRegionBED_(x) (GRegionBED(_p.opts.at(OPT_RESOURCE), x).path)
    #define GVarVCF_(x)    (GVarVCF(_p.opts.at(OPT_RESOURCE), x).path)
    #define GDecoyVCF_()   (GDecoyVCF(_p.opts.at(OPT_RESOURCE)).path)

    auto initAR = [&](UserReference &r)
    {
        const auto a1 = readAttrib(GFeatBED_(HumanAssembly::hseq));
        const auto a2 = readAttrib(GFeatBED_(HumanAssembly::hg38));
        
        r.a1 = std::shared_ptr<AttributeBed>(new AttributeBed(a1));
        r.a2 = std::shared_ptr<AttributeBed>(new AttributeBed(a2));
    };

    const auto rpath = _p.opts.at(OPT_RESOURCE);
    
    switch (_p.tool)
    {
        case Tool::RNA:
        case Tool::Meta:
        case Tool::Split:
        case Tool::Somatic:
        case Tool::Germline:
        case Tool::Calibrate:
        {
            switch (_p.tool)
            {
                case Tool::Calibrate:
                {
                    if ((!_p.opts.count(OPT_COMBINE) && !_p.opts.count(OPT_SEQUIN)) ||
                        ( _p.opts.count(OPT_COMBINE) &&  _p.opts.count(OPT_SEQUIN)))
                    {
                        throw InvalidUsageException();
                    }
                    
                    // Translation for information code
                    r.t1 = readTrans(Reader(GInfoCode(_p.opts.at(OPT_RESOURCE)).path));

                    const auto file  = _p.opts.count(OPT_SEQUIN) ? _p.opts.at(OPT_SEQUIN) : _p.opts.at(OPT_COMBINE);
                    const auto edge  = _p.opts.count(OPT_EDGE)   ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;
                    const auto build = _p.opts.count(OPT_BUILD)  ? parseBuild(_p.opts.at(OPT_BUILD)) : Detector::fromBAM(file);
                    
                    if (_p.opts.count(OPT_COMBINE))
                    {
                        const auto hr = _p.opts.count(OPT_R_HUMAN) ? _p.opts[OPT_R_HUMAN] : GRegionBED_(build);
                        const auto dr = _p.opts.count(OPT_R_DECOY) ? _p.opts[OPT_R_DECOY] : GRegionBED_(HumanAssembly::chrQS);
                        const auto rr = _p.opts.count(OPT_R_REGS)  ? _p.opts[OPT_R_REGS]  : dr;

                        if (dr == rr)
                        {
                            RegionOptions o;
                            r.r1 = readR(hr, o); // No trimming
                            r.r2 = readR(dr, o); // No trimming
                            
                            o.edge = edge;
                            r.r3 = readR(hr, o); // Trimmed
                            r.r4 = readR(dr, o); // Trimmed
                        }
                        else
                        {
                            RegionOptions o;
                            r.r1 = readR(hr, o); // No trimming
                            r.r2 = readR(dr, o); // No trimming
                            
                            o.edge = edge;
                            r.r3 = readR(hr, o); // Trimmed
                            r.r4 = readR(BedTools::intersect(dr, rr, edge), o); // Trimmed
                        }
                    }
                    else if (_p.opts.count(OPT_SEQUIN))
                    {
                        const auto hr = _p.opts.count(OPT_R_HUMAN) ? _p.opts[OPT_R_HUMAN] : GRegionBED_(build);
                        const auto rr = _p.opts.count(OPT_R_REGS)  ? _p.opts[OPT_R_REGS]  : hr;
                        
                        if (hr == rr)
                        {
                            RegionOptions o1, o2;
                            o2.edge = edge;
                            r.r1 = readR(hr, o1); // No trimming
                            r.r2 = readR(hr, o2); // Trimmed
                        }
                        else
                        {
                            RegionOptions o;
                            r.r1 = readR(hr, o); // No trimming
                            r.r2 = readR(BedTools::intersect(hr, rr, edge), o); // Trimmed
                        }
                    }

                    // Attribute regions
                    initAR(r);
                    
                    r.l1 = readL(std::bind(&Standard::addAF, &s, std::placeholders::_1), OPT_R_VCF, r, GVarVCF_(HumanAssembly::hg38));
                    r.l2 = readTSV(Reader(option(OPT_ATTTSV, GAttrBED_())), r);
                    r.l3 = readTSV(Reader(option(OPT_SYNC, GSynTSV_())), r, 2);

                    r.v1 = readV(OPT_R_VCF,  r, nullptr, GVarVCF_(HumanAssembly::hg38)); // All variants
                    r.v2 = readGV(OPT_R_VCF, r, nullptr, GVarVCF_(HumanAssembly::hg38)); // Germline variants
                    r.v3 = readSV(OPT_R_VCF, r, nullptr, GVarVCF_(HumanAssembly::hg38)); // Somatic variants
                    r.v4 = readV(OPT_R_VCF,  r, nullptr, GDecoyVCF_());                  // All variants

                    break;
                }

                case Tool::Somatic:
                case Tool::Germline:
                {
                    if ((!_p.opts.count(OPT_COMBINE) && !_p.opts.count(OPT_SEQUIN)) ||
                        ( _p.opts.count(OPT_COMBINE) &&  _p.opts.count(OPT_SEQUIN)))
                    {
                        throw InvalidUsageException();
                    }
                    
                    const auto file = _p.opts.count(OPT_SEQUIN) ? _p.opts.at(OPT_SEQUIN) : _p.opts.at(OPT_COMBINE);
                    __ha__ = _p.opts.count(OPT_BUILD) ? parseBuild(_p.opts.at(OPT_BUILD)) : Detector::fromVCF(file);

                    const auto hr   = _p.opts.count(OPT_R_HUMAN)  ? _p.opts[OPT_R_HUMAN] : GRegionBED_(__ha__);
                    const auto dr   = !_p.opts.count(OPT_COMBINE) ? hr : (_p.opts.count(OPT_R_DECOY) ? _p.opts[OPT_R_DECOY] : GRegionBED_(HumanAssembly::chrQS));
                    const auto rr   = _p.opts.count(OPT_R_REGS)   ? _p.opts[OPT_R_REGS]  : dr;
                    const auto edge = _p.opts.count(OPT_EDGE)     ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;

                    /*
                     * For combined,     a1 == chrQS and a2 == hg38
                     * For non-combined, a1 == hg38  and a2 == null
                     */
                    
                    try
                    {
                        r.a1 = std::shared_ptr<AttributeBed>(new AttributeBed(readAttrib(
                                        BedTools::intersect2(BedTools::intersect(dr, rr, 0),
                                              _p.opts.count(OPT_COMBINE) ? GFeatBED_(HumanAssembly::chrQS) : GFeatBED_(__ha__)))));
                    }
                    catch (const ZeroIntersectionError &)
                    {
                        r.a1 = std::shared_ptr<AttributeBed>(new AttributeBed());
                    }
                    
                    if (_p.opts.count(OPT_COMBINE))
                    {
                        try
                        {
                            r.a2 = std::shared_ptr<AttributeBed>(new AttributeBed(
                                        readAttrib(BedTools::intersect2(hr, GFeatBED_(__ha__)))));
                        }
                        catch (const ZeroIntersectionError &)
                        {
                            r.a2 = std::shared_ptr<AttributeBed>(new AttributeBed());
                        }
                    }

                    /*
                     * r1 = human regions without edge
                     * r2 = human regions with edge
                     * r3 = decoy regions without edge
                     * r4 = decoy regions with edge
                     */
                    
                    RegionOptions o1, o2;
                    o2.edge = edge;

                    if (dr == rr)
                    {
                        r.r1 = readR(hr, o1); // Human regions without edge
                        r.r2 = readR(hr, o2); // Human regions with edge
                        r.r3 = readR(dr, o1); // Decoy regions without edge
                        r.r4 = readR(dr, o2); // Decoy regions with edge
                    }
                    else
                    {
                        const auto f1 = BedTools::intersect(dr, rr, 0);
                        const auto f2 = BedTools::intersect(dr, rr, edge);
                        
                        if (_p.opts.count(OPT_COMBINE))
                        {
                            r.r1 = readR(hr, o1); // Human regions without edge
                            r.r2 = readR(hr, o2); // Human regions with edge
                        }
                        else
                        {
                            r.r1 = readR(f1, o1); // Human regions without edge
                            r.r2 = readR(f2, o1); // Human regions with edge
                        }
                        
                        r.r3 = readR(f1, o1); // Decoy regions without edge
                        r.r4 = readR(f2, o1); // Decoy regions with edge
                    }
                    
                    r.r5 = readR(dr, o2); // Decoy regions without intersection with edge

                    const auto vcf = _p.opts.count(OPT_COMBINE) ? GVarVCF_(HumanAssembly::chrQS) : GVarVCF_(__ha__);
                    
                    r.v1 = readV (OPT_R_VCF, r, r.r4, vcf); // All variants
                    r.v2 = readGV(OPT_R_VCF, r, r.r4, vcf); // Germline variants
                    r.v3 = readSV(OPT_R_VCF, r, r.r4, vcf); // Somatic variants

                    break;
                }

                case Tool::Meta:
                {
                    #define MMix() (MetaMix(_p.opts.at(OPT_RESOURCE)).path)

                    r.l1 = readL(std::bind(&Standard::readMMix, &s, std::placeholders::_1), OPT_MMIX, r, MMix());
                    r.l3 = readTSV(Reader(option(OPT_SYNC, GSynTSV_())), r, 2); // Unit
                    r.t1 = readTrans(Reader(GInfoCode(_p.opts.at(OPT_RESOURCE)).path));

                    break;
                }

                case Tool::RNA:
                {
                    #define RMix() (RNAMix(_p.opts.at(OPT_RESOURCE)).path)
                    
                    r.l1 = readL(std::bind(&Standard::readRMix,  &s, std::placeholders::_1), OPT_MMIX, r, RMix());
                    r.l2 = readL(std::bind(&Standard::readRGMix, &s, std::placeholders::_1), OPT_MMIX, r, RMix());
                    r.l3 = readL(std::bind(&Standard::readRLen,  &s, std::placeholders::_1), OPT_MMIX, r, RMix());
                    r.t1 = readTrans(Reader(GInfoCode(_p.opts.at(OPT_RESOURCE)).path));

                    break;
                }

                case Tool::Split:
                {
                    RegionOptions o2;
                    o2.edge = DEFAULT_EDGE;
                    
                    const auto build = _p.opts.count(OPT_BUILD) ? parseBuild(_p.opts.at(OPT_BUILD)) : HumanAssembly::hg38;
                    
                    // Attribute regions
                    initAR(r);
                    
                    r.l1 = readL(std::bind(&Standard::addAF, &s, std::placeholders::_1), OPT_R_VCF, r, GVarVCF_(build));
                    r.l2 = readTSV(Reader(option(OPT_ATTTSV, GAttrBED_())), r);
                    r.l3 = readTSV(Reader(option(OPT_SYNC, GSynTSV_())), r, 2); // Unit
                    r.t1 = readTrans(Reader(GInfoCode(_p.opts.at(OPT_RESOURCE)).path));

                    r.v1 = readV (OPT_R_VCF, r, r.r2, GVarVCF_(build)); // All variants
                    r.v2 = readGV(OPT_R_VCF, r, r.r2, GVarVCF_(build)); // Germline variants
                    r.v3 = readSV(OPT_R_VCF, r, r.r2, GVarVCF_(build)); // Somatic variants

                    /*
                     * l1 is always the base quantitative ladder. We should check if our mixture
                     * is compatible or not.
                     */
                    
                    if (_p.mix == Mix_2 && !r.l1->hasMix2())
                    {
                        throw InvalidOptionException("Mixture B specified, but only a single mixture found in reference file: " + r.l1->src);
                    }                    
                    else if (_p.mix == Mix_3 && !r.l1->hasMix3())
                    {
                        throw InvalidOptionException("Mixture C specified, but only a single mixture found in reference file: " + r.l1->src);
                    }

                    break;
                }

                default: { break; }
            }

            Standard::instance().gen.finalize(r);
            Standard::instance().rna.finalize(r);
            Standard::instance().meta.finalize(r);

            auto initSOptions = [&](SOptions &o, const FileName &ind)
            {
                assert(o.flip && o.skipMerge);

                // How many k-mers to skip?
                o.skipKM = _p.opts.count(OPT_SKIP) ? stoi(_p.opts[OPT_SKIP]) : 5;
                
                // Number of threads
                o.thr = _p.opts.count(OPT_THREAD) ? stoi(_p.opts[OPT_THREAD]) : 1;
                
                // K-mer length
                o.k = _p.opts.count(OPT_KMER) ? stoi(_p.opts[OPT_KMER]) : K_DEFAULT_L;
                
                if (_p.opts.count(OPT_NO_FLIP))
                {
                    o.flip = false;
                }
                
                o.skipMerge = !_p.opts.count(OPT_MERGE);

                // Classification rule
                o.rule = _p.opts.count(K_DEFAULT_R) ? stoi(_p.opts[K_DEFAULT_R]) : K_DEFAULT_R;
                
                o.index   = _p.opts.count(OPT_FASTA) ? _p.opts[OPT_FASTA] : ind;
                o.report  = _p.opts.count(OPT_REPORT);
                o.bam     = _p.opts.count(OPT_COMBINE);
                o.forceFQ = _p.opts.count(OPT_FASTQ);
            };
            
            auto initCalib = [&](int fKey, int aKey, const std::string &fStr, const std::string &aStr, double &calib)
            {
                if (_p.od.count(fKey))
                {
                    // Calibrate by percentage
                    calib = _p.od.at(fKey);
                    
                    if (calib < 0.0 || calib > 1.0)
                    {
                        throw std::runtime_error("Calibration fraction must be within [0,1]");
                    }
                }
                else if (_p.od.count(aKey))
                {
                    // Calibrate by absolute
                    if ((calib = _p.od.at(aKey)) < 1)
                    {
                        throw std::runtime_error("Number of reads must not be less than 1");
                    }
                }
                else
                {
                    // If no calibration specified, it must be -1
                    assert(calib == NO_CALIBRATION);
                }
            };
            
            switch (_p.tool)
            {
                case Tool::Meta:
                {
                    MSplit::Options o;
                    initSOptions(o, MetaFA(_p.opts.at(OPT_RESOURCE)).path);

                    o.mix = _p.mix;
                    assert(o.oneS == NO_CALIBRATION && o.oneL == NO_CALIBRATION && o.secL == NO_CALIBRATION);

                    double tmp = -1;
                    
                    initCalib(OPT_S_FRAC, OPT_S_ABSO, "--sequin_fraction", "--sequin_absolute", o.oneS);
                    initCalib(OPT_L_FRAC, OPT_L_ABSO, "--ladder_fraction", "--ladder_absolute", tmp);
                    
                    if (tmp > 1) { o.oneL = tmp; }
                    else         { o.secL = tmp; }
                    
                    start<MSplit>("meta", [&](const MSplit::Options &)
                    {
                        if (o.bam) { MSplit::report(_p.opts[OPT_COMBINE], "", o); }
                        else       { MSplit::report(_p.opts[OPT_1], _p.opts.count(OPT_2) ? _p.opts[OPT_2] : _p.opts[OPT_1], o); }
                    }, o);

                    break;
                }
                    
                case Tool::RNA:
                {
                    RSplit::Options o;
                    initSOptions(o, RNAFA(_p.opts.at(OPT_RESOURCE)).path);
                    
                    o.mix = _p.mix;
                    assert(o.oneS == NO_CALIBRATION && o.oneL == NO_CALIBRATION && o.secL == NO_CALIBRATION);
                    
                    initCalib(OPT_S_FRAC, OPT_S_ABSO, "--sequin_fraction", "--sequin_absolute", o.oneS);
                    
                    start<RSplit>("rna", [&](const RSplit::Options &)
                    {
                        if (o.bam) { RSplit::report(_p.opts[OPT_COMBINE], "", o); }
                        else       { RSplit::report(_p.opts[OPT_1], _p.opts.count(OPT_2) ? _p.opts[OPT_2] : _p.opts[OPT_1], o); }
                    }, o);

                    break;
                }
                    
                case Tool::Split:
                {
                    GSplit::Options o;
                    initSOptions(o, GSeqFA_());
                    
                    if (o.bam) { o.forceFQ = true; }

                    start<GSplit>("split", [&](const GSplit::Options &)
                    {
                        if (o.bam) { GSplit::report(_p.opts[OPT_COMBINE], "", o); }
                        else       { GSplit::report(_p.opts[OPT_1], _p.opts.count(OPT_2) ? _p.opts[OPT_2] : _p.opts[OPT_1], o); }
                    }, o);

                    break;
                }
                    
                case Tool::Somatic:
                case Tool::Germline:
                {
                    GVariant::Options o;
                    
                    o.build    = __ha__;
                    o.combined = _p.opts.count(OPT_COMBINE);
                    o.edge     = _p.opts.count(OPT_EDGE) ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;
                    o.uBED     = _p.opts.count(OPT_R_REGS) ? _p.opts[OPT_R_REGS] : MISSING;
                    o.report   = _p.opts.count(OPT_REPORT);
                    
                    const auto f1 = o.combined ? _p.opts.at(OPT_COMBINE) : _p.opts.count(OPT_SAMPLE) ? _p.opts.at(OPT_SAMPLE) : "";
                    const auto f2 = o.combined ? "" : _p.opts.at(OPT_SEQUIN);
                    
                    if (_p.tool == Tool::Germline)
                    {
                        start<GGerm>(o.base = "germline", [&](const GVariant::Options &)
                        {
                            GGerm::report(f1, f2, o);
                        }, o);
                    }
                    else
                    {
                        start<GSomatic>(o.base = "somatic", [&](const GVariant::Options &)
                        {
                            GSomatic::report(f1, f2, o);
                        }, o);
                    }

                    break;
                }

                case Tool::Calibrate:
                {
                    typedef GCalibrate::Method Method;
                    
                    GCalibrate::Options o;
                    initSOptions(o, GSeqFA_());
                    o.edge = _p.opts.count(OPT_EDGE) ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;                    
                    
                    o.decoy = GSeqDecoy(_p.opts.at(OPT_RESOURCE)).path;
                    o.nonDecoy = GSeqFA_();
                    
                    // How to calibrate?
                    const auto meth = _p.opts.count(OPT_METHOD) ? _p.opts[OPT_METHOD] : "mean";
                    
                    if      (meth == "mean")   { o.meth = Method::Mean;   }
                    else if (meth == "median") { o.meth = Method::Median; }
                    else if (isFloat(meth))    { o.meth = Method::Percent; o.p = stof(meth); }
                    else
                    {
                        throw std::runtime_error("Unknown method: " + meth);
                    }

                    if (o.bam)
                    {
                        o.cMode = GCalibrate::CalibrateMode::Combined;
                        analyze_1<GCalibrate>("calibrate", OPT_COMBINE, o);
                    }
                    else
                    {
                        if (o.meth != Method::Percent && !_p.opts.count(OPT_SAMPLE))
                        {
                            throw MissingOptionError("--sample");
                        }
                        
                        o.cMode = GCalibrate::CalibrateMode::TwoBAM;
                        analyze_2<GCalibrate>("calibrate", o.meth != Method::Percent ? OPT_SAMPLE : OPT_SEQUIN, OPT_SEQUIN, o);
                    }
                    
                    break;
                }

                default : { break; }
            }

            break;
        }

        default : { break; }
    }
}

extern int parse_options(int argc, char ** argv)
{
    char cwd[1024];
    
    auto printError = [&](const std::string &x)
    {
        std::cerr << "***********************" << std::endl;
        std::cerr << "[ERRO]: " << x << std::endl;
        std::cerr << "***********************" << std::endl << std::endl;
    };
    
    if (getcwd(cwd, sizeof(cwd)))
    {
        __working__ = cwd;
    }
    
    try
    {
        parse(argc, argv);
        return 0;
    }
    catch (const FailedCommandException &ex)
    {
        printError(std::string(ex.what()));
    }
    catch (const InvalidFormatException &ex)
    {
        printError("Invalid file format: " + std::string(ex.what()));
    }
    catch (const InvalidUsageException &)
    {
        printError("Invalid usage. Please check and try again.");
    }
    catch (const InvalidToolError &ex)
    {
        printError("Invalid command. Unknown tool: " + ex.v + ". Please check your usage and try again.");
    }
    catch (const InvalidOptionException &ex)
    {
        printError((boost::format("Invalid usage. Unknown option: %1%") % ex.o).str());
    }
    catch (const InvalidValueException &ex)
    {
        printError((boost::format("Invalid command. %1% not expected for %2%.") % ex.v % ex.o).str());
    }
    catch (const MissingOptionError &ex)
    {
        const auto format = "Invalid command. Mandatory option is missing. Please specify %1%.";
        printError((boost::format(format) % ex.o).str());
    }
    catch (const MissingBundleError &ex)
    {
        printError((boost::format("%1%%2%") % "Failed to find bundle at: " % ex.path).str());
    }
    catch (const InvalidBundleError &ex)
    {
        printError((boost::format("%1%%2%") % "Invalid bundle. Missing bundle folder: " % ex.file).str());
    }
    catch (const InvalidFileError &ex)
    {
        printError((boost::format("%1%%2%") % "Invalid command. File is invalid: " % ex.file).str());
    }
    catch (const InvalidDependency &ex)
    {
        printError(ex.what());
    }
    catch (const std::runtime_error &ex)
    {
        printError(ex.what());
    }

    return 1;
}

int main(int argc, char ** argv)
{
    srand(time(0));
    return parse_options(argc, argv);
}
