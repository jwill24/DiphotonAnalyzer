#ifndef DiphotonAnalyzer_EventAnalyzer_RPCalibrationConstants_h
#define DiphotonAnalyzer_EventAnalyzer_RPCalibrationConstants_h

#include "DataFormats/Provenance/interface/RunID.h"

namespace CTPPSAlCa
{
    typedef struct RPCalibrationConstants {
        RPCalibrationConstants() :
            x_disp_l_n( 0. ), x_disp_l_f( 0. ),
            x_disp_r_n( 0. ), x_disp_r_f( 0. ) {}
        RPCalibrationConstants( float l_n, float l_f, float r_n, float r_f ) :
            x_disp_l_n( l_n ), x_disp_l_f( l_f ),
            x_disp_r_n( r_n ), x_disp_r_f( r_f ) {}

        // dispersions (x-axis)
        float x_disp_l_n, x_disp_l_f;
        float x_disp_r_n, x_disp_r_f;
    } RPCalibrationConstants;

    RPCalibrationConstants getCalibrationConstants( const edm::RunNumber_t& );
}

#endif
