#include "data/reader.hpp"
#include "tools/attributes.hpp"

using namespace Anaquin;

AttributeBed Anaquin::readAttrib(const FileName &f)
{
    AttributeBed x;
    
    ParserBed::parse(Reader(x.src = f), [&](const ParserBed::Data &d, Progress)
    {
        std::vector<Token> toks;
        split(d.name, "_", toks);
        
        if (!x.count(toks[0]))
        {
            x[toks[0]] = std::shared_ptr<BedData>(new BedData());
        }

        x.vals[toks[0]].insert(remove(d.name, toks[0] + "_"));
        (*x[toks[0]])[d.cID].l2d[d.l] = d;
    });
    
    return x;
}
