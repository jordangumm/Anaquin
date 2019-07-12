#ifndef FILE_WRITER_HPP
#define FILE_WRITER_HPP

#include <iomanip>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <boost/format.hpp>
#include "writers/writer.hpp"

namespace Anaquin
{
    class FileWriter : public Writer<>
    {
        public:

            FileWriter(const Path &path) : path(path) {}
            ~FileWriter() { close(); }
        
            static void create(const Path &path, const FileName &file, const std::string &txt)
            {
                FileWriter w(path);
                w.open(file);
                w.write(txt);
                w.close();
            }

            inline void close() override
            {
                if (_o)
                {
                    _o->close();
                    _o.reset();
                    _o = nullptr;
                }
            }

            inline void open(const FileName &file) override
            {
                isScript = isEnd(file, ".R") || isEnd(file, ".py");

                if (!path.empty())
                {
                    createD(path);
                }
                
                const auto target = !path.empty() ? path + "/" + file : file;
                _o = std::shared_ptr<std::ofstream>(new std::ofstream(target));
                
                if (!_o->good())
                {
                    throw std::runtime_error("Failed to open: " + target);
                }
            }

            inline void write(const std::string &x) override
            {
                if (isScript)
                {
                    *(_o) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << trim(x);
                }
                else
                {
                    *(_o) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << x;
                }
                
                *(_o) << std::endl;
            }

        private:
        
            // Trim script?
            bool isScript;
        
            std::string path;
            std::shared_ptr<std::ofstream> _o;
    };
}

#endif
