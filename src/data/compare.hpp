#ifndef COMPARE_HPP
#define COMPARE_HPP

#include <math.h>

namespace Anaquin
{
    /*
     * This class represents a data-wrapper for Cuffcompare
     */

    struct Compare
    {
        // Metrics at the base level
        double b_sp, b_sn;

        // Metrics at the exon level
        double e_sp, e_sn, e_fsp, e_fsn;
        
        // Metrics at the intron level
        double i_sp = NAN, i_sn = NAN, i_fsp = NAN, i_fsn = NAN;
        
        // Metrics at the intron-chain level
        double c_sp = NAN, c_sn = NAN, c_fsp = NAN, c_fsn = NAN;

        // Metrics at the locus level
        double l_sp, l_sn, l_fsp, l_fsn;

        // Metrics at the transcript level
        double t_sp, t_sn, t_fsp, t_fsn;

        double novelExonsP, novelIntronsP;
        double missedExonsP, missedIntronsP;

        unsigned novelExonsN,  novelExonsR,    novelIntronsN, novelIntronsR;
        unsigned missedExonsR, missedIntronsR, missedExonsN,  missedIntronsN;
    };
}

#endif