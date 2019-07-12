#include "tools/tools.hpp"
#include "data/reference.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

struct IntersectResults
{
    std::set<SequinID> diffs, inters;
};

template <typename T> IntersectResults intersect(const std::set<T> &t1, const std::set<T> &t2)
{
    std::set<SequinID> x, y;
    
    for (const auto &i : t1) { x.insert(static_cast<SequinID>(i)); }
    for (const auto &i : t2) { y.insert(static_cast<SequinID>(i)); }
    
    A_ASSERT(!x.empty() && !y.empty());
    
    IntersectResults c;
    
    std::set_intersection(x.begin(),
                          x.end(),
                          y.begin(),
                          y.end(),
                          std::inserter(c.inters, c.inters.begin()));

    std::set_difference(x.begin(),
                        x.end(),
                        y.begin(),
                        y.end(),
                        std::inserter(c.diffs, c.diffs.begin()));

    std::set_difference(y.begin(),
                        y.end(),
                        x.begin(),
                        x.end(),
                        std::inserter(c.diffs, c.diffs.begin()));

    return c;
}

/*
 * ------------------------- Genomics Analysis -------------------------
 */

struct GenomicsRef::GenomicsRefImpl
{
    // Empty Implementation
};

GenomicsRef::GenomicsRef() : _impl(new GenomicsRefImpl()) {}

Counts GenomicsRef::nType(std::shared_ptr<VCFLadder> v, Variation x) const
{
    return v->data.count_(x);
}

void GenomicsRef::validate(const UserReference &r)
{
    _t1 = r.t1;
    _l1 = r.l1; _l2 = r.l2; _l3 = r.l3;
    _r1 = r.r1; _r2 = r.r2; _r3 = r.r3; _r4 = r.r4; _r5 = r.r5;
    _v1 = r.v1; _v2 = r.v2; _v3 = r.v3;
    _a1 = r.a1; _a2 = r.a2;
}
