#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include "data/locus.hpp"

namespace Anaquin
{
    struct Alignment
    {
        operator const Locus &() const { return l; }

        // Primary alignment (not std::string for effiency)
        ChrID cID;

        // Location of the alignment
        Locus l;
        
        // If this field is false, no assumption can be made to other fields
        bool mapped;
        
        // Mapping quality
        int mapq;

        // Bitwise FLAG
        int flag;

        /*
         * SAM flag fields
         */

        bool isReverseC;        
        bool isPaired;
        bool isAllAligned;
        bool isAligned;
        bool isMateAligned;
        bool isMateReverse;
        bool isFirstPair;
        bool isSecondPair;
        bool isDuplicate;
        bool isPrimary;
        bool isSupplement;        
        bool isPassed;
        
        // Secondary alignment? Typically used for alternative mappings when multiple mappings are presented
        bool isSecondary;
        
        // Signed observed template length
        Base tlen;
        
        /*
         * Optional fields (for efficiency they are not reset by the parser)
         */
        
        // Eg: B7_591:6:155:12:674
        ReadName name;

        // Segment sequence
        Sequence seq;
        
        // ASCII of base QUALity plus 33
        std::string qual;        
    };
}

#endif
