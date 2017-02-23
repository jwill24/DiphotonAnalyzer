#ifndef DiphotonAnalyzer_EventAnalyzer_ProtonKinematicsUtils_h
#define DiphotonAnalyzer_EventAnalyzer_ProtonKinematicsUtils_h

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "DiphotonAnalyzer/EventAnalyzer/interface/RPAlignmentConstants.h"

#include <utility> // std::pair

namespace ProtonUtils
{
  typedef std::pair<unsigned int, const TotemRPLocalTrack&> localtrack_t;

  float tracksDistance( const CTPPSAlCa::RPAlignmentConstants&, const localtrack_t& near_p, const localtrack_t& far_p );
}

#endif
