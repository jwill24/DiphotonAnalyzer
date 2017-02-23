#include "DiphotonAnalyzer/EventAnalyzer/interface/ProtonKinematicsUtils.h"

#include <iostream>

namespace ProtonUtils
{
  float
  tracksDistance( const CTPPSAlCa::RPAlignmentConstants& align, const localtrack_t& near_p, const localtrack_t& far_p )
  {
    CTPPSAlCa::RPAlignmentConstants::Quantities align_near = align.quantities( near_p.first ),
                                                align_far = align.quantities( far_p.first );
    /*std::cout << "---> " << align_near.x << "/" << align_near.y << std::endl
              << "     " << near_p.second.getX0() << "/" << near_p.second.getY0() << std::endl;*/
    const float z_far = far_p.second.getZ0();
    const TVector2 far_xy_ext = near_p.second.getTrackPoint( z_far ) + TVector2( align_near.x, align_near.y ),
                   far_xy_obs = TVector2( far_p.second.getX0(), far_p.second.getY0() ) + TVector2( align_far.x, align_far.y );
    return ( far_xy_ext-far_xy_obs ).Mod() / 10.; // mm -> cm
  }
}
