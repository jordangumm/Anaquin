#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include "data/standard.hpp"
#include "writers/writer.hpp"

namespace Anaquin
{
    struct LimitStats
    {
        // Absolute detection limit
        Limit limit;
    };

    class WriterOptions
    {
        public:
            WriterOptions() : writer(std::make_shared<MockWriter>()),
                              logger(std::make_shared<MockWriter>()),
                              output(std::make_shared<MockWriter>()) {}
        
            WriterOptions(const WriterOptions &x)
            {
                this->work = x.work;
                this->writer = x.writer;
                this->logger = x.logger;
                this->output = x.output;
            }
        
            // Working directory
            Path work;

            std::shared_ptr<Writer<>> writer, logger, output;

            inline void wait(const std::string &s) const
            {
                log("[WAIT]: " + s);
                out("[WAIT]: " + s);
            }

            inline void logWait(const std::string &s) const
            {
                log("[WAIT]: " + s);
            }
              
            inline void warn(const std::string &s) const
            {
                log("[WARN]: " + s);
                out("[WARN]: " + s);
            }

            inline void logWarn(const std::string &s) const
            {
                log("[WARN]: " + s);
            }

            inline void analyze(const std::string &s) const
            {
                info("Analyzing: " + s);
            }
        
            inline void generate(const FileName &f) const
            {
                info("Generating " + f);
            }
        
            inline void info(const std::string &s) const
            {
                log("[INFO]: " + s);
                out("[INFO]: " + s);
            }
        
            inline void logInfo(const std::string &s) const
            {
                log("[INFO]: " + s);
            }

            inline void error(const std::string &s) const
            {
                info("[ERROR]: " + s);
            }
        
        private:
        
            inline void out(const std::string &s) const { if (output) { output->write(s); } }
            inline void log(const std::string &s) const { if (logger) { logger->write(s); } }
    };

    template <typename Build> struct AnalyzerOptions : public WriterOptions
    {
        AnalyzerOptions() {}
        AnalyzerOptions(const AnalyzerOptions &x) : WriterOptions(x), name(x.name), cmd(x.cmd), build(x.build) {}
        
        // Eg: "split"_summary.stats
        std::string name;
        
        // Full comamnd
        std::string cmd;

        // Human hg38, gencode etc
        Build build;        
    };
}

#endif
