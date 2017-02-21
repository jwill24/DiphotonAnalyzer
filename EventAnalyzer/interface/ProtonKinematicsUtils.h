#ifndef DiphotonAnalyzer_EventAnalyzer_ProtonKinematicsUtils_h
#define DiphotonAnalyzer_EventAnalyzer_ProtonKinematicsUtils_h

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

namespace ProtonUtils
{
    float tracksDistance( const TotemRPLocalTrack& near_p, const TotemRPLocalTrack& far_p );
}

#endif
