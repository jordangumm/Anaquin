#ifndef REFERENCE_HPP
#define REFERENCE_HPP

#include "data/vData.hpp"
#include "data/bData.hpp"
#include "data/ladder.hpp"
#include "data/reader.hpp"
#include "data/variant.hpp"
#include "data/dinters.hpp"
#include "tools/attributes.hpp"

namespace Anaquin
{
    struct VCFLadder
    {
        FileName src;
        
        VCFData data;
        
        // Allele frequency ladder
        Ladder lad;
    };

    enum class Tool
    {
        RNA,
        Meta,
        Split,
        Calibrate,
        Germline,
        Somatic
    };
    
    struct Ladder;
    struct VCFLadder;

    typedef std::map<std::string, std::string> Translation;
    
    struct UserReference
    {
        std::shared_ptr<AttributeBed> a1, a2;
        
        std::shared_ptr<Ladder> l1, l2, l3, l4;

        std::shared_ptr<VCFLadder> v1, v2, v3;
        
        std::shared_ptr<Translation> t1;
        
        std::shared_ptr<BedData> r1, r2, r3, r4, r5;
    };
    
    class Reference
    {
        public:

            typedef std::set<SequinID> SequinIDs;

            inline Concent input1(const SequinID &x, Mixture m = Mix_1) const { return _l1->input(x, m); }
            inline Concent input2(const SequinID &x, Mixture m = Mix_1) const { return _l2->input(x, m); }

            // Allele frequency from VCF reference (not from mixture file)
            inline Proportion af(const SequinID &x) const { return _v1->lad.input(x, Mix_1); }

            inline std::shared_ptr<AttributeBed> a1() const { return _a1; }
            inline std::shared_ptr<AttributeBed> a2() const { return _a2; }
        
            inline std::shared_ptr<Ladder> l1() const { return _l1; }
            inline std::shared_ptr<Ladder> l2() const { return _l2; }
            inline std::shared_ptr<Ladder> l3() const { return _l3; }

            inline std::shared_ptr<BedData> r1() const { return _r1; }
            inline std::shared_ptr<BedData> r2() const { return _r2; }
            inline std::shared_ptr<BedData> r3() const { return _r3; }
            inline std::shared_ptr<BedData> r4() const { return _r4; }
            inline std::shared_ptr<BedData> r5() const { return _r5; }

            inline std::shared_ptr<VCFLadder> v1() const { return _v1; }
            inline std::shared_ptr<VCFLadder> v2() const { return _v2; }
            inline std::shared_ptr<VCFLadder> v3() const { return _v3; }

            inline std::shared_ptr<Translation> t1() const { return _t1; }
        
            inline void finalize(const UserReference &r)
            {
                validate(r);
            }

        protected:

            virtual void validate(const UserReference &) = 0;

            std::shared_ptr<Translation> _t1;
            std::shared_ptr<AttributeBed> _a1, _a2;
        
            // Sequin regions
            std::shared_ptr<BedData> _r1, _r2, _r3, _r4, _r5;

            // VCF references and allele frequency ladder
            std::shared_ptr<VCFLadder> _v1, _v2, _v3;
        
            // Sequin ladders
            std::shared_ptr<Ladder> _l1, _l2, _l3;
    };

    class GenomicsRef : public Reference
    {
        public:
            GenomicsRef();

            // Number of references by variation
            Counts nType(std::shared_ptr<VCFLadder>, Variation) const;

        protected:
            void validate(const UserReference &) override;

        private:
            struct GenomicsRefImpl;
            std::shared_ptr<GenomicsRefImpl> _impl;
    };
    
    class RnaRef : public Reference
    {
        protected:
            void validate(const UserReference &r) override { _l1 = r.l1; _l2 = r.l2; _l3 = r.l3; _t1 = r.t1; }
    };
    
    class MetaRef : public Reference
    {
        protected:
            void validate(const UserReference &r) override { _l1 = r.l1; _l2 = r.l2; _l3 = r.l3; _t1 = r.t1; }
    };}

#endif
