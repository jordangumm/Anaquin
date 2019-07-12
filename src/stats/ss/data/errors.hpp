#ifndef SS_ERRORS_HPP
#define SS_ERRORS_HPP

#include <stdexcept>

namespace SS
{
    #define SS_ASSERT(cond, message) \
        if (!(cond)) { \
            throw std::runtime_error(message); \
        }
}

#endif