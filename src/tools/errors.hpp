#ifndef ERRORS_HPP
#define ERRORS_HPP

#include <string>
#include <assert.h>
#include <stdexcept>

namespace Anaquin
{
    #define A_ASSERT(cond) \
    assert(cond);

    #define A_CHECK(cond, message) \
    if (!(cond)) { throw std::runtime_error(message); }

    #define A_THROW(message) \
    throw std::runtime_error(message);
    
    typedef std::exception ZeroIntersectionError;

    struct InvalidFileError : public std::exception
    {
        InvalidFileError(const std::string &file) : file(file) {}
        const std::string file;
    };

    struct FailedCommandException : public std::runtime_error
    {
        FailedCommandException(const std::string &msg) : std::runtime_error(msg) {}
    };
    
    struct InvalidFormatException : public std::runtime_error
    {
        InvalidFormatException(const std::string &msg) : std::runtime_error(msg) {}
    };
}

#endif
