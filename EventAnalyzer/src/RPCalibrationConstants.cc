#include "DiphotonAnalyzer/EventAnalyzer/interface/RPCalibrationConstants.h"

namespace CTPPSAlCa
{
    RPCalibrationConstants
    getCalibrationConstants( const edm::RunNumber_t& )
    {
        RPCalibrationConstants cc;
        cc.x_disp_l_f = 9.22e-2;
        cc.x_disp_l_n = 9.26e-2;
        cc.x_disp_r_f = 5.16e-2;
        cc.x_disp_r_n = 5.81e-2;
        return cc;
    }
}
