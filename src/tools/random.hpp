#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <map>
#include <string>
#include <memory>
#include <klib/khash.h>
#include "tools/errors.hpp"

namespace Anaquin
{
    class RandomSelection
    {
        public:

            RandomSelection(double prob, int seed = rand())  : _seed(seed), _prob(prob)
            {
                A_ASSERT(prob >= 0.0);
            }

            inline bool select(const std::string &hash) const
            {
                const uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(hash.c_str()) ^ _seed);
                return ((double)(k&0xffffff) / 0x1000000 >= _prob);
            }

        private:

            // Random seed
            int _seed;

            // The probability of selection
            const double _prob;
    };

    typedef std::map<std::string, std::shared_ptr<RandomSelection>> Selection;
}

#endif
