#ifndef ATTRIBUTES_HPP
#define ATTRIBUTES_HPP

#include <set>
#include <map>
#include <sstream>
#include "data/bData.hpp"
#include "tools/tools.hpp"

namespace Anaquin
{
    typedef std::string AttribKey;
    typedef std::string AttribGrp;
    typedef std::string AttribStr;

    inline AttribKey attrKey(const AttribStr &x) { return first(x, "_"); }
    inline AttribGrp attrGrp(const AttribStr &x) { return last(x, "_");  }
    
    class AttributeBed : public std::map<AttribKey, std::shared_ptr<BedData>>
    {
        public:

            inline const std::set<AttribKey> &keys() const
            {
                if (_keys.empty())
                {
                    for (const auto &i : (*this))
                    {
                        _keys.insert(i.first);
                    }
                }
            
                return _keys;
            }

            inline std::string strForKeys() const
            {
                std::stringstream ss;
                for (const auto &i : (*this)) { ss << "\t" << i.first; }            
                return ss.str();
            }

            inline std::string strForNulls() const
            {
                if (_nulls.empty())
                {
                    std::stringstream ss;
                    for (auto i = 0u; i < keys().size(); i++) { ss << "\t"; ss << MISSING; }
                    _nulls = ss.str();
                }
                
                return _nulls;
            }

            static std::string strForVals(const std::map<AttribKey, AttribGrp> &am)
            {
                std::stringstream ss;
                for (const auto &i : am) { ss << "\t"; ss << i.second; }
                return ss.str();
            }
        
            inline std::string strForLocus(const ChrID &cID, const Locus &l)
            {
                std::stringstream ss;
                
                for (const auto &i : keys())
                {
                    ss << "\t"; ss << strForLocus(i, cID, l);
                }
                
                return ss.str();
            }
        
            inline std::string strForLocus(const AttribKey &key, const ChrID &cID, const Locus &l)
            {
                if (!_cache.count(key))
                {
                    _cache[key] = (*this)[key]->ginters();
                }
                
                const auto &m = _cache[key];
                std::vector<GInterval *> o;
            
                if (m.count(cID) && m.at(cID).overlap(l, &o))
                {
                    std::vector<std::string> os;
                
                    std::transform(o.begin(), o.end(), std::back_inserter(os), [&](GInterval *i)
                    {
                        return noFirst(i->name(), "_");
                    });

                    return os.back();
                }
            
                return MISSING;
            }
        
            FileName src;
            std::map<AttribKey, std::set<AttribGrp>> vals;

        private:
            mutable std::string _nulls;
            mutable std::set<AttribKey> _keys;
            mutable std::map<AttribKey, std::map<ChrID, GIntervals<>>> _cache;
    };
    
    AttributeBed readAttrib(const FileName &);
}

#endif
