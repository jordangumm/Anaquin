#ifndef SS_INTERNAL_TEST_HPP
#define SS_INTERNAL_TEST_HPP

#include <functional>
#include <ss/stats.hpp>
#include <ss/data/errors.hpp>
#include <ss/data/results.hpp>
#include <ss/internal/dist.hpp>

namespace SS
{
    namespace Internal
    {
        template <typename Iter> static Results tTestOneSample
                (const Iter &x, Real h0, P conf = 0.95, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            SS_ASSERT(x.size(), "Samples size must not be zero");
            
            const auto u  = mean(x);
            const auto n  = x.size();
            const auto s  = getSD(x) / sqrt(n);
            const auto t  = (u - h0) / s;
            const auto df = n - 1;
            
            const auto c = Internal::critical(std::bind(qt, _1, df), conf, type);
            const auto p = Internal::pval(t, pt(t, df), type);
            
            const auto lc = (type != Less    ? u - (s * fabs(c)) : -INFINITY);
            const auto uc = (type != Greater ? u + (s * fabs(c)) :  INFINITY);
            
            Results r;
            
            r.addStats(t);
            r.addP(p);
            r.addLCI(lc);
            r.addUCI(uc);
            r.addDF(df);
            
            return r;
        }
        
        template <typename Iter> static Results tTestTwoSamples
                (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.95, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            SS_ASSERT(x.size() == y.size(), "Samples must match in dimension");
            SS_ASSERT(x.size(), "Samples size must not be zero");
            
            const auto u1 = mean(x);
            const auto u2 = mean(y);
            const auto s1 = getSD(x);
            const auto s2 = getSD(y);
            const auto n1 = x.size();
            const auto n2 = y.size();
            const auto u  = u1 - u2;
            
            const DF df = 2.0 * n1 - 2.0;
            
            const auto s  = sqrt((((n1 - 1.0) * s1 * s1) + ((n2 - 1.0) * s2 * s2)) / (n1 + n2 - 2.0));
            const auto se = s * sqrt(1.0/n1 + 1.0/n2);
            
            const auto t = (u1 - u2 - h0) / se;
            const auto c = Internal::critical(std::bind(qt, _1, df), conf, type);
            const auto p = Internal::pval(t, pt(t, df), type);
            
            const auto lc = type != Less    ? (u1 - u2) - fabs(c) * se : -INFINITY;
            const auto uc = type != Greater ? (u1 - u2) + fabs(c) * se :  INFINITY;
            
            Results r;
            
            r.addStats(t);
            r.addP(p);
            r.addLCI(lc);
            r.addUCI(uc);
            r.addDF(df);
            
            return r;
        }
        
        template <typename Iter> static Results tTestPaired
        (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.05, TestType type = TwoSided)
        {
            SS_ASSERT(x.size() == y.size(), "Samples must match in dimension");
            SS_ASSERT(x.size(), "Samples size must not be zero");
            
            std::vector<Real> diffs;
            diffs.resize(x.size());
            
            for (auto i = 0; i < x.size(); i++)
            {
                diffs[i] = x[i] - y[i];
            }
            
            return tTestOneSample(diffs, h0, conf, type);
        }
        
        template <typename Iter> static Results tTestWelch
        (const Iter &x, const Iter &y, Real h0 = 0.0, P conf = 0.05, TestType type = TwoSided)
        {
            using namespace std::placeholders;
            
            SS_ASSERT(x.size() == y.size(), "Samples must match in dimension");
            SS_ASSERT(x.size(), "Samples size must not be zero");
            
            const auto u1 = mean(x);
            const auto u2 = mean(y);
            const auto s1 = getSD(x);
            const auto s2 = getSD(y);
            const auto n1 = x.size();
            const auto n2 = y.size();
            
            const auto s = std::sqrt((s1 * s1 / n1) + (s2 * s2 / n2));
            const auto t = (u1 - u2 - h0) / s;
            
            const DF df = std::pow(s1 * s1 / n1 + s2 * s2 / n2, 2) /
            (std::pow(s1 * s1 / n1, 2) / (n1 - 1) +
             std::pow(s2 * s2 / n2, 2) / (n2 - 1));
            
            const auto c = Internal::critical(std::bind(qt, _1, df), conf, type);
            const auto p = Internal::pval(t, pt(t, df), type);
            
            const auto lc = type != Less    ? (u1 - u2) - fabs(c) * s : -INFINITY;
            const auto uc = type != Greater ? (u1 - u2) + fabs(c) * s :  INFINITY;
            
            Results r;
            
            r.addStats(t);
            r.addP(p);
            r.addLCI(lc);
            r.addUCI(uc);
            r.addDF(df);
            
            return r;
        }
    }
}

#endif