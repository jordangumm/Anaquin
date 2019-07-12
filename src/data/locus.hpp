#ifndef LOCUS_HPP
#define LOCUS_HPP

#include <string>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <stdexcept>
#include "data/data.hpp"

namespace Anaquin
{
    class Locus
    {
        public:

            Locus(const Locus &l1, const Locus &l2, const std::string &key = "") : _key(key)
            {
                end   = std::max(l1.end,   l2.end);
                start = std::min(l1.start, l2.start);
            }

            Locus(Base start = 0, Base end = 0, const std::string &key = "")
                    : start(start), end(end), _key(key)
            {
                if (end < start)
                {
                    throw std::runtime_error("Locus: end < start");
                }
            }

            inline std::string key() const
            {
                return !_key.empty() ? _key : std::to_string(start) + "-" + std::to_string(end);
            }
        
            template <typename T, typename R, template <typename, typename = std::allocator<T>> class Inputs>
                    static std::vector<R> merge(const Inputs<T> &x)
            {
                std::vector<R> sorted;

                // std::copy doesn't allow implicit type conversion...
                for (const auto &i : x)
                {
                    sorted.push_back(i);
                }

                std::sort(sorted.begin(), sorted.end(), [&](const Locus &l1, const Locus &l2)
                {
                    return (l1.start < l2.start) || (l1.start == l2.start && l1.end < l2.end);
                });
                
                std::vector<R> merged;
                
                for (auto i = 0; i < sorted.size();)
                {
                    // We'll need the index for skipping i
                    auto j = i + 1;

                    auto super = sorted[i];
                    
                    // Look forward until the end of an overlap region
                    for (; j < sorted.size(); j++)
                    {
                        if (static_cast<Locus>(sorted[j]).overlap(super))
                        {
                            super += sorted[j];
                        }
                        else
                        {
                            break;
                        }
                    }

                    // Construct the super-loci for the region
                    merged.push_back(super);

                    i = j;
                }

                return merged;
            }

            void merge(const Locus &l)
            {
                end   = std::max(end, l.end);
                start = std::min(start, l.start);
            }

            void inter(const Locus &l)
            {
                end   = std::min(end, l.end);
                start = std::max(start, l.start);
            }
                                                
            inline void add(const std::string &key) { _key = key; }

            inline Base length() const { return (end - start + 1); }

            inline Base overlap(const Locus &l) const
            {
                if (l.start > end || start > l.end)
                {
                    return 0;
                }
                else if (start <= l.start && end >= l.end)
                {
                    return l.length();
                }
                else if (start >= l.start && end <= l.end)
                {
                    return length();
                }
                else if (end >= l.end)
                {
                    return l.end - start + 1;
                }
                else
                {
                    return end - l.start + 1;
                }
            }

            inline bool contains(const Locus &q) const
            {
                return (q.start >= start && q.end <= end);
            }

            inline void operator+=(const Locus &l)
            {
                this->start = std::min(start, l.start);
                this->end   = std::max(end, l.end);
            }
            
            inline void operator+=(long n) { start += n; end += n; }
            inline void operator-=(long n) { start -= n; end -= n; }

            inline bool operator!=(const Locus &l) const { return !operator==(l); }
            inline bool operator==(const Locus &l) const { return start == l.start && end == l.end; }

            inline Locus operator+(const Locus &l) const
            {
                return Locus(std::min(start, l.start), std::max(end, l.end));
            }

            inline bool operator<(const Locus &l) const
            {
                return start < l.start || (start == l.start && end < l.end);
            }

            inline operator std::string() const
            {
                return std::to_string(start) + "-" + std::to_string(end);
            }

            Base start, end;

        private:
            std::string _key;
    };
}

#endif
