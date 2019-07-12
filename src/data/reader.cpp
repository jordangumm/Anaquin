#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>

#include <zlib.h>
#include <boost/algorithm/string.hpp>

#include "data/reader.hpp"

namespace Anaquin
{
    class GzipFileReader : public AbstractReader
    {
    public:
        GzipFileReader() = delete;
        GzipFileReader(const GzipFileReader&) = delete;
        
        GzipFileReader(const std::string& fname);
        ~GzipFileReader();
        
        virtual bool nextLine(std::string& s) const;
        virtual void reset();
        
    private:
        static const size_t BUFLEN = 1024 * 1024;
        
        gzFile_s *data = nullptr;
        std::shared_ptr<char> buf;
    };

    static const std::string GZIPSuffix = ".gz";
    
    bool
    AbstractReader::nextTokens(std::vector<std::string>& tokens, const std::string& c) const
    {
        tmp.clear();
        
        if (nextLine(tmp)) {
            tokens.clear();
            boost::split(tokens, tmp, boost::is_any_of(c));
            std::for_each(std::begin(tokens), std::end(tokens), [this](std::string& s) { trim(s); });
            return true;
        }
        
        return false;
    }
    
    std::string
    AbstractReader::src() const
    {
        return fname;
    }
    
    void
    AbstractReader::trim(std::string& s) const
    {
        boost::trim(s);
    }
    
    TxtFileReader::TxtFileReader(const std::string& fname)
    : AbstractReader(fname)
    {
        data = std::make_shared<std::ifstream>(fname, std::ios::in);
        if (data->fail()) {
            throw std::runtime_error("couldn't open file");
        }
    }
    
    bool
    TxtFileReader::nextLine(std::string& s) const
    {
        auto status = !!std::getline(*data, s);
        trim(s);
        
        return status;
    }
    
    void
    TxtFileReader::reset()
    {
        data->clear();
        data->seekg(0, std::ios::beg);
    }
    
    GzipFileReader::GzipFileReader(const std::string& fname)
    : AbstractReader(fname)
    {
        data = gzopen(fname.c_str(), "rb");
        if (!data) {
            throw std::runtime_error("couldn't open file");
        }
        
        buf = std::shared_ptr<char>(new char[BUFLEN], [](char* p) { delete[] p; });
    };
    
    GzipFileReader::~GzipFileReader()
    {
        if (data) {
            auto rc = gzclose(data);
            if (rc != Z_OK) {
                // TODO
            }
        }
    }
    
    bool
    GzipFileReader::nextLine(std::string& s) const
    {
        auto p = gzgets(data, buf.get(), BUFLEN);
        if (!p) {
            return false;
        }
        
        s.assign(p);
        trim(s);
        
        return true;
    }
    
    void
    GzipFileReader::reset()
    {
        gzrewind(data);
    }
    
    Reader::Reader(const std::string& fname, bool forceGZ)
    {
        if (!fname.empty()) {
            if (fname.length() >= GZIPSuffix.length()
                && fname.compare(fname.length() - GZIPSuffix.length(), GZIPSuffix.length(), GZIPSuffix) == 0) {
                type = FileType::GZIP;
            }
            
            if (forceGZ) {
                type = FileType::GZIP;
            }
            
            if (type == FileType::GZIP) {
                reader = std::make_shared<GzipFileReader>(fname);
            } else {
                reader = std::make_shared<TxtFileReader>(fname);
            }
        }
    };
    
} // namespace
