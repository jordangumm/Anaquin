#ifndef BUNDLE_HPP
#define BUNDLE_HPP

#include <map>
#include <glob.h>
#include <vector>
#include "data/data.hpp"
#include "tools/tools.hpp"

namespace Anaquin
{
    struct Resource
    {
        FileName path;
    };

    struct Bundle
    {
        static std::vector<std::string> myGlob(const std::string &p)
        {
            using namespace std;
            vector<string> names;

            glob_t gr;
            memset(&gr, 0, sizeof(gr));
            
            int r = glob(p.c_str(), GLOB_TILDE, NULL, &gr);
            if (r != 0)
            {
                globfree(&gr);
                return names;
            }
            
            for (std::size_t i = 0; i < gr.gl_pathc; ++i)
            {
                names.push_back(string(gr.gl_pathv[i]));
            }
            
            globfree(&gr);
            return names;
        }
        
        static FileName latest(const FileName &s1, const FileName &s2)
        {
            std::map<float, FileName> m;
            
            for (auto file : myGlob(s1 + "*" + s2))
            {
                // Eg: "2.5"
                auto v = remove(remove(remove(remove(file, s1), s2), "beta"), "_");
                
                // Don't crash for unexpected files
                try { m[stof(v)] = file; } catch (...) {}
            }
            
            const auto x = m.empty() ? "" : m.rbegin()->second;
            
            if (x.empty() || !exists(x))
            {
                throw std::runtime_error("Failed to find resource bundle for " + (s1 + "*" + s2));
            }

            return x;
        }
    };
}

#endif
