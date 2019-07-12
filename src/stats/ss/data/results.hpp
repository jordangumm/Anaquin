#ifndef SS_RESULTS_HPP
#define SS_RESULTS_HPP

#include <map>
#include <string>
#include <ss/data/data.hpp>

namespace SS
{
    namespace Keys
    {
        typedef std::string Key;
        
        /*
         * Keys for all statistical tests
         */
        
        const Key FStats = "fstats";
        const Key PValue = "pvalue";
        
        /*
         * Keys related to ANOVA
         */
        
        const Key ErrorDF = "errorDF";
        const Key ErrorSS = "errorSS";
        const Key ErrorMS = "errorMS";

        const Key ModelDF = "modelDF";
        const Key ModelSS = "modelSS";
        const Key ModelMS = "modelMS";

        const Key TotalDF = "totalDF";
        const Key TotalSS = "totalSS";
        const Key TotalMS = "totalMS";
    }

    const auto KeyDF    = "DF";
    const auto KeyP     = "PV";
    const auto KeyStats = "Stats";
    const auto KeyLCI   = "LowerCI";
    const auto KeyUCI   = "UpperCI";
    
    struct Results : public std::map<std::string, Real>
    {
        inline void addP(Real x)     { (*this)[KeyP]     = x; }
        inline void addDF(Real x)    { (*this)[KeyDF]    = x; }
        inline void addStats(Real x) { (*this)[KeyStats] = x; }
        
        inline void addLCI(Real x) { (*this)[KeyLCI] = x; }
        inline void addUCI(Real x) { (*this)[KeyUCI] = x; }

        // P-value
        inline Real p() const { return (*this).at(KeyP); }
        
        // Degree of freedom
        inline Real df() const { return (*this).at(KeyDF); }

        // Lower confidence interval
        inline Real lci() const { return (*this).at(KeyLCI); }
        
        // Upper confidence interval
        inline Real uci() const { return (*this).at(KeyUCI); }

        // Test statistics
        inline Real stats() const { return (*this).at(KeyStats); }
    };

    class StatsResults : public std::map<std::string, Real>
    {
        public:
        
            typedef std::string StatsName;

            /*
             * Create an object for keeping statistical results by string keys.
             */

            StatsResults(const StatsName &name) : _name(name) {}

        private:

            // Name of the statistics (eg: "lm", "anova")
            StatsName _name;
    };
}

#endif