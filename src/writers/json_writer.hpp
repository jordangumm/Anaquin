#ifndef JSON_WRITER_HPP
#define JSON_WRITER_HPP

#include <map>
#include <string>
#include <sstream>
#include "writers/file_writer.hpp"
#include <boost/algorithm/string/replace.hpp>

namespace Anaquin
{
    struct JSONWriter : public FileWriter
    {
        JSONWriter(const Path &path) : FileWriter(path) {}
        
        inline void write(const std::map<std::string, std::string> &x)
        {
            std::stringstream ss;
            ss << "{ ";
            
            for (const auto &i : x)
            {
                ss << "\"" << i.first << "\":\"" << i.second << "\",";
            }
            
            ss << "}";
            
            auto json = ss.str();
            boost::replace_all(json, ",}", "}");
            FileWriter::write(json);
        }
    };
}

#endif
