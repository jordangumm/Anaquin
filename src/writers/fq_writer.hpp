#ifndef FQ_WRITER_HPP
#define FQ_WRITER_HPP

#include <fstream>
#include "tools/tools.hpp"

namespace Anaquin
{
    class FQWriter
    {
        public:
            FQWriter(const FileName &file)
            {
                _o = std::shared_ptr<std::ofstream>(new std::ofstream(file, std::ios::binary | std::ios::out));
            }

            ~FQWriter() { close(); }
        
            void write(const ReadName &r, const Sequence &s, const Sequence &q)
            {
                std::stringstream str;
                str << "@" + std::string(r) + "\n";
                str << s + "\n";
                str << "+\n";
                str << std::string(q) + "\n";
                (*_o) << compressGZ(str.str());
            };
        
            void close()
            {
                if (_o)
                {
                    _o->close();
                    _o = nullptr;
                }
            }

        private:
            std::shared_ptr<std::ofstream> _o;
    };
}

#endif
