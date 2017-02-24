#include "DiphotonAnalyzer/EventAnalyzer/interface/ProtonKinematicsUtils.h"

#include <iostream>

namespace ProtonUtils
{
  float
  tracksDistance( const CTPPSAlCa::RPAlignmentConstants& align, const localtrack_t& near_p, const localtrack_t& far_p )
  {
    if ( near_p.first/100!=far_p.first/100 ) { // avoid dealing with opposite arms
      return -1.;
    }

    // start by extracting the RP alignment constants for the corresponding pots
    const CTPPSAlCa::RPAlignmentConstants::Quantities align_near = align.quantities( near_p.first ),
                                                      align_far = align.quantities( far_p.first );

    const TVector2 al_near( align_near.x, align_near.y ),
                   al_far( align_far.x, align_far.y );

    const TotemRPLocalTrack near_track( near_p.second ),
                            far_track( far_p.second);

    const TVector2 far_xy_ext = near_track.getTrackPoint( far_track.getZ0() ) + al_near,
                   far_xy_obs = TVector2( far_track.getX0(), far_track.getY0() ) + al_far;
    return ( far_xy_ext-far_xy_obs ).Mod() / 10.; // mm -> cm
  }
}
