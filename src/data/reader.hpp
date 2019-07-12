#ifndef READER_HPP
#define READER_HPP

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include <unistd.h>

namespace Anaquin
{
    
    class AbstractReader
    {
    public:
        AbstractReader() = delete;
        AbstractReader(const AbstractReader& r) = delete;
        
        AbstractReader(const std::string& fname)
        : fname(fname){};
        
        ~AbstractReader()
        {
        }
        
        virtual void reset() = 0;
        
        // Returns the next line in the file
        virtual bool nextLine(std::string&) const = 0;
        
        // Returns the next line and parse it into tokens
        bool nextTokens(std::vector<std::string>& tokens, const std::string& c) const;
        
        // Returns description for the source
        std::string src() const;
        
    protected:
        std::string fname;
        
        void trim(std::string& s) const;
        
    private:
        mutable std::string tmp;
    };
    
    class TxtFileReader : public AbstractReader
    {
    public:
        TxtFileReader() = delete;
        TxtFileReader(const TxtFileReader&) = delete;
        
        TxtFileReader(const std::string& fname);
        
        ~TxtFileReader()
        {
        }
        
        // Returns the next line in the file
        virtual bool nextLine(std::string& s) const;
        virtual void reset();
        
    private:
        std::shared_ptr<std::ifstream> data;
    };
    
    class Reader
    {
    public:
        Reader() = delete;
        
        Reader(const std::string &, bool forceGZ = false);
        ~Reader(){};
        
        Reader(const Reader &r)
        {
            reader = r.reader;
            reader->reset();
        }
        
        static bool valid(const std::string &file)
        {
            std::ifstream r(file);
            return r.good() && r.peek() != std::ifstream::traits_type::eof();
        }
        
        void reset()
        {
            reader->reset();
        }
        
        std::string src() const
        {
            return reader->src();
        }
        
        // Returns the next line in the file
        bool nextLine(std::string& s) const
        {
            return reader->nextLine(s);
        };
        
        // // Returns the next line and parse it into tokens
        bool nextTokens(std::vector<std::string>& tokens, const std::string& c) const
        {
            return reader->nextTokens(tokens, c);
        };
        
    private:
        enum class FileType {
            TXT,
            GZIP,
        };
        
        std::shared_ptr<AbstractReader> reader;
        FileType type = FileType::TXT;
    };
} // namespace

#endif
