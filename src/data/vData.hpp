#ifndef VCF_DATA_HPP
#define VCF_DATA_HPP

#include "tools/tools.hpp"
#include "parsers/parser_vcf.hpp"

namespace Anaquin
{
    struct VData
    {
        std::map<Base, Variant> b2v;
        std::map<Variation, std::set<Variant>> m2v;
    };

    class VCFData : public std::map<ChrID, VData>
    {
        public:
        
            inline const Variant *findByKey(long key)
            {
                lazyInitKey2Var();
                return _key2Var->count(key) ? (*_key2Var)[key] : nullptr;
            }
        
            inline const Variant *find(const SequinID &x)
            {
                for (const auto &i : *this)
                {
                    for (const auto &j : i.second.b2v)
                    {
                        if (j.second.name == x)
                        {
                            return &j.second;
                        }
                    }
                }
            
                return nullptr;
            }
        
            inline std::map<SequinID, const Variant *> seq2Var() const
            {
                std::map<SequinID, const Variant *> x;
            
                for (const auto &i : *this)
                {
                    for (const auto &j : i.second.m2v)
                    {
                        for (const auto &k : j.second)
                        {
                            x[k.name] = &k;
                        }
                    }
                }
            
                return x;
            }

            inline std::set<Variant> vars() const
            {
                std::set<Variant> x;
            
                for (const auto &i : *this)
                {
                    for (const auto &j : i.second.m2v)
                    {
                        for (const auto &k : j.second)
                        {
                            x.insert(k);
                        }
                    }
                }
            
                return x;
            }

            inline const Variant * findVar(const ChrID &x, const Locus &l)
            {
                if (!count(x))
                {
                    return nullptr;
                }
                else if (at(x).b2v.count(l.start))
                {
                    return &(at(x).b2v.at(l.start));
                }

                return nullptr;
            }

            inline Counts count_(const ChrID &x, Variation m) const
            {
                return count(x) && at(x).m2v.count(m) ? at(x).m2v.at(m).size() : 0;
            }

            inline Counts count_(Variation m) const
            {
                return countMap(*this, [&](const ChrID &x, const VData &)
                {
                    return count_(x, m);
                });
            }
        
        private:

            void lazyInitKey2Var()
            {
                if (!_key2Var)
                {
                    _key2Var = std::shared_ptr<std::map<long, const Variant *>>(new std::map<long, const Variant *>());
                    
                    for (const auto &i : *this)
                    {
                        for (const auto &j : i.second.b2v)
                        {
                            (*_key2Var)[j.second.key()] = &(j.second);
                        }
                    }
                }
            }
        
            std::shared_ptr<std::map<long, const Variant *>> _key2Var;
    };
}

#endif
