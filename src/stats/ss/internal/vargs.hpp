#ifndef SS_INTERNAL_VARGS_HPP
#define SS_INTERNAL_VARGS_HPP

#include <vector>

namespace SS
{
    namespace Internal
    {
        template <typename T, typename F> void vArgs(F f, std::vector<const T *> &p, const T &t)
        {
            p.push_back(&t);
            f(p);
        }

        template<typename T, typename F, typename... Args> void vArgs(F f, std::vector<const T *> &p, const T &t, Args... args)
        {
            p.push_back(&t);
            vArgs(f, p, args...) ;
        }

        template<typename T, typename F, typename... Args> void vArgs(F f, const T &t, Args... args)
        {
            std::vector<const T *> p;
            vArgs(f, p, t, args...) ;
        }
    }
}

#endif