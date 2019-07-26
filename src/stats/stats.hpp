#ifndef STATS_HPP
#define STATS_HPP

#include <set>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "data/data.hpp"
#include "stats/linear.hpp"

namespace Anaquin
{
    template <typename T1, typename T2> typename T1::value_type quant(const T1 &x, T2 q)
    {
        assert(q >= 0.0 && q <= 1.0);
        
        auto y = std::vector<typename T1::value_type>(x.begin(), x.end());
        auto const i = q * x.size();
        
        std::nth_element(y.begin(), y.begin() + i, y.end());
        return (y.size() % 2 || y.size() == 1) ? y.at(i) : 0.5 * (y.at(i-1) + y.at(i));
    }

    template <typename T1, typename T2> typename T1::value_type hackQuant(const T1 &x, T2 q)
    {
        assert(q >= 0.0 && q <= 1.0);
        auto const i = q * x.size();
        return x.at(i);
    }

    template <typename T> typename T::value_type med(const T &x)
    {
        return quant(x, 0.5);
    }

    template <typename T> typename T::value_type min(const T &x)
    {
        return *(std::min_element(x.begin(), x.end()));
    }
    
    template <typename T> typename T::value_type max(const T &x)
    {
        return *(std::max_element(x.begin(), x.end()));
    }
    
    Counts RFilterR(const FileName &, const FileName &, const std::set<Label> &);
    
    // Filter columns
    void RFilterC(const FileName &, const FileName &, const std::set<Label> &, bool keep = false);
    
    typedef std::function<double (const std::vector<double> &)> Apply;
    std::map<double, double> RBinaryTSV(const FileName &, const Label &, const Label &);
    
    enum class Imputation
    {
        None,
        ToZero,
        Remove
    };
    
    void RAggregate(const FileName &, const FileName &, const Label &, Apply, Imputation impute = Imputation::None);
    void RAggregateSD(const FileName &, const FileName &, const Label &, Imputation impute = Imputation::None);
    void RAggregateCount(const FileName &, const FileName &, const Label &);
    void RAggregateMean(const FileName &, const FileName &, const Label &, Imputation impute = Imputation::None);
    void RAggregateSum(const FileName &, const FileName &, const Label &, Imputation impute = Imputation::None);
    void RAggregateMax(const FileName &, const FileName &, const Label &);
    void RAggregateMin(const FileName &, const FileName &, const Label &);
    void RAggregateQ0(const FileName &, const FileName &, const Label &);
    void RAggregateQ25(const FileName &, const FileName &, const Label &);
    void RAggregateQ50(const FileName &, const FileName &, const Label &);
    void RAggregateQ75(const FileName &, const FileName &, const Label &);
    void RAggregateQ100(const FileName &, const FileName &, const Label &);

    typedef std::function<std::string (const std::string &)> RApplyF;
    void RApply(const FileName &, const FileName &, const Label &, RApplyF);

    // Do we keep it when there's a match?
    void RGrep(const FileName &, const FileName &, const Label &, const std::string &, bool keep = true, bool isNum = false);

    double RSum(const FileName &, const Label &);
    Counts RCount(const FileName &, const Label &, const std::string &);

    bool RHead(const FileName &, const Label &);
    
    void RMeanCV(const FileName &, const FileName &, const Label &, const Label &, Counts nExp = 2, Counts Obs = 2, Counts nCV = 2);

    // Returns mean ratios for a ladder table
    double RLadTable(const FileName &, const FileName &, const Label &);

    SequinStats RLinear(const FileName &src, const Label &, const Label &, const Label &);
}

#endif
