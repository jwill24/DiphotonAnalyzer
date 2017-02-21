#include "DiphotonAnalyzer/EventAnalyzer/interface/ProtonKinematicsUtils.h"

namespace ProtonUtils
{
    float
    tracksDistance( const TotemRPLocalTrack& near_p, const TotemRPLocalTrack& far_p )
    {
        const float z_far = far_p.getZ0();
        const TVector2 far_xy_ext = near_p.getTrackPoint( z_far ), far_xy_obs = TVector2( far_p.getX0(), far_p.getY0() );
        return (far_xy_ext-far_xy_obs).Mod()/1.e3; // mm -> m
    }
}
