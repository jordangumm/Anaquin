#ifndef __SAMPLE_HPP__
#define __SAMPLE_HPP__

#include <iomanip>
#include <sstream>
#include <boost/format.hpp>
#include "data/data.hpp"
#include "stats/ss/stats.hpp"
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

namespace Anaquin
{
    static inline std::string d2str(double x)
    {
        std::ostringstream out;
        out << std::setprecision(6) << x;
        return out.str();
    }
    
    template <typename T> static std::string concat(const std::vector<T> &x, std::string (*f)(T) = d2str)
    {
        return boost::algorithm::join(x | boost::adaptors::transformed(static_cast<std::string(*)(T)>(f)), ", ");
    }
    
    static std::string concat(const std::vector<std::string> &x, const std::string &join = "\'")
    {
        return (join + boost::algorithm::join(x, join + ", " + join) + join);
    }

    #define STRING(x) static_cast<std::string>(x)

    template <typename T> class SSamples
    {
        public:
        
            inline void add(const T &x)
            {
                _data.push_back(x);
            }
        
            inline std::size_t size() const { return _data.size(); }
        
            virtual operator std::string() const = 0;

        protected:
        
            std::vector<T> _data;
    };
    
    struct SReals : public SSamples<double>
    {
        inline double mean() const
        {
            auto data = SSamples<double>::_data;
            
            // Some samples might not define...
            data.erase(std::remove_if(data.begin(), data.end(),[](double x) { return std::isnan(x); }), data.end());
            
            return SS::mean(data);
        }
        
        virtual operator std::string() const
        {
            auto data = SSamples<double>::_data;
            
            // Some samples might not define...
            data.erase(std::remove_if(data.begin(), data.end(),[](double x) { return std::isnan(x); }), data.end());
            
            if (data.size() == 0)
            {
                return "-";
            }
            else if (data.size() > 1)
            {
                return (boost::format("%1$.2f \u00B1 %2$.2f") % SS::mean(data) % SS::SD(data)).str();
            }
            else
            {
                return (boost::format("%1$.2f") % data.front()).str();
            }
        }
    };
    
    struct SStrings : public SSamples<std::string>
    {
        virtual operator std::string() const
        {
            return concat(SSamples<std::string>::_data, "");
        }
    };

    struct SCounts : public SSamples<Counts>
    {
        virtual operator std::string() const
        {
            const auto &data = SSamples<Counts>::_data;
            
            if (data.empty())
            {
                return "-";
            }
            else if (data.size() > 1)
            {
                return (boost::format("%1% \u00B1 %2%") % SS::mean(data) % SS::SD(data)).str();
            }
            else
            {
                return (boost::format("%1%") % data.front()).str();
            }
        }
    };
}

#endif
