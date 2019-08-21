#include <rapidcsv.h>
#include "tools/tools.hpp"
#include "stats/stats.hpp"
#include "stats/ss/stats.hpp"

using namespace Anaquin;

template <typename T> bool isMissing(const T& x) { return x == MISSING || x == "." || x == "-"; }

static std::vector<Token> RGetColumns(rapidcsv::Document &doc)
{
    std::vector<Token> toks;
    
    for (auto i = 0; i < doc.GetColumnCount(); i++)
    {
        toks.push_back(doc.GetColumnName(i));
    }
    
    return toks;
}

static void RWriteHeader(std::ofstream &w, const std::vector<Token> &toks)
{
    for (auto i = 0; i < toks.size(); i++)
    {
        w << toks.at(i) << (i != toks.size() - 1 ? "\t" : "");
    }
}

Counts Anaquin::RFilterR(const FileName &src, const FileName &dst, const std::set<std::string> &cols)
{
    std::ofstream w(dst);
    
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    RWriteHeader(w, heads);
    
    std::vector<Token> tmp;
    tmp.resize(heads.size());
    
    std::map<std::string, std::size_t> m;
    
    for (auto iter = cols.begin(); iter != cols.end(); iter++)
    {
        // Guarantee there's always an entry
        m[*iter] = std::find(heads.begin(), heads.end(), *iter) - heads.begin();
    }
    
    // Number of filtered
    auto k = 0;
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        
        auto shouldWrite = [&]()
        {
            for (auto iter = m.begin(); iter != m.end(); iter++)
            {
                if (isMissing(row[iter->second]))
                {
                    return false;
                }
            }
            
            return true;
        };
        
        if (shouldWrite())
        {
            k++;
            w << std::endl;
            
            for (auto j = 0; j < row.size(); j++)
            {
                w << row[j] << (j != row.size() - 1 ? "\t" : "");
            }
        }
    }
    
    w.close();
    return k;
}

SequinStats Anaquin::RLinear(const FileName &src, const Label &name, const Label &x, const Label &y)
{
    SequinStats s;
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    const auto i = std::find(heads.begin(), heads.end(), x) - heads.begin();    // Index for independent
    const auto j = std::find(heads.begin(), heads.end(), y) - heads.begin();    // Index for dependent
    const auto n = std::find(heads.begin(), heads.end(), name) - heads.begin(); // Index for ID

    assert(i != heads.size());
    assert(j != heads.size());
    assert(n != heads.size());

    for (auto k = 0; k < doc.GetRowCount(); k++)
    {
        const auto row = doc.GetRow<std::string>(k);
        
        if (!isMissing(row[i]) && !isMissing(row[j]))
        {
            s.add(row[n], stod(row[i]), stod(row[j]));
        }
    }
    
    return s;
}

std::map<double, double> Anaquin::RBinaryTSV(const FileName &file, const Label &key, const Label &val)
{
    std::map<double, double> r;
    
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto cols = RGetColumns(doc);
    assert(cols.size() == 2);
    
    auto ki = std::find(cols.begin(), cols.end(), key) - cols.begin();
    auto vi = std::find(cols.begin(), cols.end(), val) - cols.begin();

    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        assert(row.size() == 2);
        r[stod(row[ki])] = (row[vi] == MISSING ? NAN : stod(row[vi]));
    }
    
    return r;
}

void Anaquin::RMeanCV(const FileName &src, const FileName &dst, const Label &key, const Label &val, Counts nExp, Counts nObs, Counts nCV)
{
    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();
    const auto tmp3 = tmpFile();
    const auto tmp4 = tmpFile();

    RFilterC(src,  tmp1, std::set<std::string> { key, val }, true);
    RFilterR(tmp1, tmp2, std::set<std::string> { key });
    RFilterR(tmp2, tmp1, std::set<std::string> { val });
    RAggregateMean(tmp1, tmp2, key);  // tmp2 holds the results
    RAggregateSD(tmp1, tmp3, key);    // tmp3 holds the results
    RAggregateCount(tmp1, tmp4, key); // tmp4 holds the results

    auto m = RBinaryTSV(tmp2, key, val); // Mean
    auto s = RBinaryTSV(tmp3, key, val); // SD
    auto n = RBinaryTSV(tmp4, key, val); // Counts
    
    assert(m.size() == s.size());
    assert(m.size() == n.size());
    
    std::ofstream w(dst);
    w << "NAME\tCOUNT\tMEAN\tCV";

    for (const auto &i : m)
    {
        assert(s.count(i.first));
        assert(n.count(i.first));
        
        // Coefficient of variation
        const auto cv = !s[i.first] ? "NA" : toString((double) s[i.first] / m[i.first], nCV);

        w << std::endl << toString(i.first, nExp) << "\t" << n.at(i.first) << "\t" << toString(m[i.first], nObs) << "\t" << cv;
    }
    
    w.close();
}

void Anaquin::RAggregate(const FileName &src, const FileName &dst, const Label &col, Apply f, Imputation impute)
{
    std::ofstream w(dst);
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    
    auto cols = RGetColumns(doc);
    assert(cols.size() >= 2);
    RWriteHeader(w, cols);
    
    // Column index
    auto k = std::find(cols.begin(), cols.end(), col) - cols.begin();
    
    std::set<Label>  k1;
    std::set<double> k2;

    std::map<Label, std::map<Label, std::vector<double>>> m;

    for (auto i = 0; i < cols.size(); i++)
    {
        if (i != k) { m[cols[i]]; }
    }

    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        assert(row.size() >= 2);
        
        for (auto j = 0; j < row.size(); j++)
        {
            if (j != k)
            {
                auto x = row[j];
                
                if (x == MISSING)
                {
                    switch (impute)
                    {
                        case Imputation::None:   { break; }
                        case Imputation::ToZero: { x = "0"; break; }
                        case Imputation::Remove: { continue; }
                    }
                }
                
                k1.insert(row[k]); if (isNumber(row[k])) { k2.insert(stod(row[k])); }
                m[cols[j]][row[k]].push_back(stod(x));
            }
        }
    }
    
    std::vector<Label> keys(k1.begin(), k1.end());
    
    if (k1.size() == k2.size())
    {
        std::sort(keys.begin(), keys.end(), [] (const std::string &x, const std::string &y) {
            return std::stod(x) < std::stod(y);
        });
    }

    for (const auto &key : keys)
    {
        w << std::endl << key;

        for (const auto &col : cols)
        {
            if (m.count(col))
            {
                if (m[col].count(key))
                {
                    w << "\t" << f(m[col][key]);
                }
                else
                {
                    w << "\t" << MISSING;
                }
            }
        }
    }

    w.close();
}

double Anaquin::RLadTable(const FileName &src, const FileName &dst, const Label &col)
{
    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();

    // Only ladders
    RGrep(src, tmp1, col, "_LD_");

    // Writer for adding "Copy" column
    std::ofstream w1(tmp2);

    rapidcsv::Document doc(tmp1, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    w1 << "COPY\t"; RWriteHeader(w1, heads);
    
    // Column index for the ID
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();

    assert(k != heads.size());
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        w1 << std::endl;

        auto cn = -1;
        if (isSubstr(row[k], "_A")) { cn = 1; }
        if (isSubstr(row[k], "_B")) { cn = 2; }
        if (isSubstr(row[k], "_C")) { cn = 4; }
        if (isSubstr(row[k], "_D")) { cn = 8; }
        assert(cn != -1);
        
        w1 << cn;
        for (auto j = 0; j < row.size(); j++)
        {
            w1 << "\t" << row[j];
        }
    }
    
    // tmp2 has the results
    w1.close();
    
    const auto key = "COPY";
    const auto val = "Q50";
    
    RFilterC(tmp2, tmp1, std::set<Label> { key, val }, true);
    const auto nonZero = !RFilterR(tmp1, tmp2, std::set<Label> { val }); // tmp2 holds the results

    std::map<double, double> mu, sd, q0, q25, q50, q75, q100;
    
    if (!nonZero)
    {
        RAggregateMean(tmp2, tmp1, key); mu   = RBinaryTSV(tmp1, key, val);
        RAggregateSD  (tmp2, tmp1, key); sd   = RBinaryTSV(tmp1, key, val);
        RAggregateMin (tmp2, tmp1, key); q0   = RBinaryTSV(tmp1, key, val);
        RAggregateQ25 (tmp2, tmp1, key); q25  = RBinaryTSV(tmp1, key, val);
        RAggregateQ50 (tmp2, tmp1, key); q50  = RBinaryTSV(tmp1, key, val);
        RAggregateQ75 (tmp2, tmp1, key); q75  = RBinaryTSV(tmp1, key, val);
        RAggregateMax (tmp2, tmp1, key); q100 = RBinaryTSV(tmp1, key, val);
    }
    
    assert(mu.size() == sd.size());
    assert(sd.size() == q0.size());
    assert(q0.size() == q25.size());
    assert(q25.size() == q50.size());
    assert(q50.size() == q75.size());
    assert(q75.size() == q100.size());

    std::ofstream w2(dst);
    w2 << "COPY\tMEAN\tSD\tCV\tQ0\tQ25\tQ50\tQ75\tQ100\tRATIO";

    std::vector<double> keys;
    for (const auto &i : mu) { keys.push_back(i.first); }

    // List of ratios
    std::vector<double> ros;
    
    for (auto i = 0; i < keys.size(); i++)
    {
        const auto cv = mu[keys[i]] == 0.0 ? "NA" : toString((double) sd[keys[i]] / mu[keys[i]], 4.0);
        const auto ro = !i ? NAN : q50[keys[i]] / q50[keys[i-1]];
        const auto cn = std::to_string((int) std::pow(2, i));
        
        #define X(x) "\t" << x[keys[i]]
        w2 << std::endl << cn << X(mu) << X(sd) << "\t" << cv << X(q0) << X(q25) << X(q50) << X(q75) << X(q100) << "\t" << replaceNA(ro, 4);
        
        if (!std::isnan(ro))
        {
            ros.push_back(ro);
        }
    }
    
    w2.close();
    
    // Mean ratios
    return ros.empty() ? NAN : SS::mean(ros);
}

void Anaquin::RApply(const FileName &src, const FileName &dst, const Label &col, RApplyF f)
{
    std::ofstream w(dst);
    
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    RWriteHeader(w, heads);
    
    // Column index
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        auto row = doc.GetRow<std::string>(i);
        row[k] = f(row[k]);
        
        w << std::endl;
        
        for (auto j = 0; j < row.size(); j++)
        {
            w << row[j] << (j != row.size() - 1 ? "\t" : "");
        }
    }
    
    w.close();
}

bool Anaquin::RHead(const FileName &file, const Label &col)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    return (std::find(heads.begin(), heads.end(), col) - heads.begin()) < heads.size();
}

Counts Anaquin::RCount(const FileName &file, const Label &col, const std::string &x)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    assert(k < heads.size());

    Counts n = 0;
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        if (doc.GetRow<std::string>(i)[k] == x)
        {
            n++;
        }
    }
    
    return n;
}

double Anaquin::RSum(const FileName &file, const Label &col)
{
    rapidcsv::Document doc(file, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    assert(k < heads.size());
    
    double n = 0;

    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto x = doc.GetRow<std::string>(i)[k];
        
        if (x != MISSING)
        {
            n += stod(x);
        }
    }
    
    return n;
}

void Anaquin::RGrep(const FileName &src, const FileName &dst, const Label &col, const std::string &val, bool keep, bool isNum)
{
    std::ofstream w(dst);
    
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    RWriteHeader(w, heads);
    
    auto k = std::find(heads.begin(), heads.end(), col) - heads.begin();
    assert(k < heads.size());
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        const auto match = (!isNum && isSubstr(row[k], val)) || (isNum && stod(row[k]) == stod(val));
        
        if ((keep && match) || (!keep && !match))
        {
            w << std::endl;
            
            for (auto j = 0; j < row.size(); j++)
            {
                w << row[j] << (j != row.size() - 1 ? "\t" : "");
            }
        }
    }
    
    w.close();
}

void Anaquin::RAggregateMean(const FileName &src, const FileName &dst, const Label &col, Imputation impute)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return SS::mean(x);
    }, impute);
}

void Anaquin::RAggregateCount(const FileName &src, const FileName &dst, const Label &col)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return std::count_if(x.begin(), x.end(), [&](double x)
        {
            return x > 0 && !std::isinf(x) && !std::isnan(x);
        });
    });
}

void Anaquin::RAggregateSD(const FileName &src, const FileName &dst, const Label &col, Imputation impute)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return SS::SD(x);
    }, impute);
}

void Anaquin::RAggregateMin(const FileName &src, const FileName &dst, const Label &col)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return min(x);
    });
}

void Anaquin::RAggregateQ25(const FileName &src, const FileName &dst, const Label &col)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return quant(x, 0.25);
    });
}

void Anaquin::RAggregateQ50(const FileName &src, const FileName &dst, const Label &col)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return med(x);
    });
}

void Anaquin::RAggregateQ75(const FileName &src, const FileName &dst, const Label &col)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return quant(x, 0.75);
    });
}

void Anaquin::RAggregateMax(const FileName &src, const FileName &dst, const Label &col)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return max(x);
    });
}

void Anaquin::RAggregateSum(const FileName &src, const FileName &dst, const Label &col, Imputation impute)
{
    return RAggregate(src, dst, col, [&](const std::vector<double> &x)
    {
        return std::accumulate(x.begin(), x.end(), 0.0);
    }, impute);
}

void Anaquin::RFilterC(const FileName &src, const FileName &dst, const std::set<Label> &cols, bool keep)
{
    rapidcsv::Document doc(src, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams('\t'));
    auto heads = RGetColumns(doc);
    
    std::set<std::size_t> m;
    
    for (auto iter = cols.begin(); iter != cols.end(); iter++)
    {
        m.insert(std::find(heads.begin(), heads.end(), *iter) - heads.begin());
    }
    
    std::ofstream w(dst);
    std::vector<Token> tmp;
    
    for (auto i = 0; i < heads.size(); i++)
    {
        if ((!keep && !m.count(i)) || (keep && m.count(i)))
        {
            tmp.push_back(heads[i]);
        }
    }
    
    for (auto i = 0; i < tmp.size(); i++)
    {
        w << tmp[i] << (i != tmp.size() - 1 ? "\t" : "");
    }
    
    for (auto i = 0; i < doc.GetRowCount(); i++)
    {
        const auto row = doc.GetRow<std::string>(i);
        
        auto started = false;
        for (auto j = 0; j < row.size(); j++)
        {
            if (!j)
            {
                w << std::endl;
            }
            
            if ((!keep && !m.count(j)) || (keep && m.count(j)))
            {
                w << (started ? "\t" : "") << row[j];
                started = true;
            }
        }
    }
    
    w.close();
}

