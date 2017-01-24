#include "DiphotonAnalyzer/EventAnalyzer/interface/ProtonKinematicsUtils.h"

namespace ProtonUtils
{
    void defaultAlignmentConstants()
    {
        // all alignment constants are given in mm
        RPAlignmentConstants period1, period2;
        period1.x_shift_l_f = -3.40;
        period1.x_shift_l_n = -0.90;
        period1.x_shift_r_f = -2.40;
        period1.x_shift_r_n = -2.75;
        period2.x_shift_l_f = -3.90;
        period2.x_shift_l_n = -1.45;
        period2.x_shift_r_f = -2.85;
        period2.x_shift_r_n = -3.25;
        alignmentMap[0] = period1;
        alignmentMap[274244] = period2;
    }
        //                                                                                                                                         }
        //                                                                                                                                                 return ac;
        //

    void
    fillAlignmentConstants( const char* filename )
    {
    }

    RPAlignmentConstants
    getAlignmentConstants( const edm::RunNumber_t& run_id )
    {
        for ( std::map<edm::RunNumber_t,RPAlignmentConstants>::const_iterator it=alignmentMap.begin()+1; it!=alignmentMap.end(); it++ ) {
            if ( ( ( it-1 ).first<run_id ) and ( it.first>=run_id ) ) return it.second;
        }
        return RPAlignmentConstants();
    }
}
