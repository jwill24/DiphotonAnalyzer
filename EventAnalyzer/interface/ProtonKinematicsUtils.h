#ifndef DiphotonAnalyzer_EventAnalyzer_ProtonKinematicsUtils_h
#define DiphotonAnalyzer_EventAnalyzer_ProtonKinematicsUtils_h

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "DataFormats/Provenance/interface/RunID.h"

namespace ProtonUtils
{
    float tracksDistance( const TotemRPLocalTrack& near_p, const TotemRPLocalTrack& far_p )
    {
        const float z_far = far_p.getZ0();
        const TVector2 far_xy_ext = near_p.getTrackPoint( z_far ), far_xy_obs = TVector2( far_p.getX0(), far_p.getY0() );
        return (far_xy_ext-far_xy_obs).Mod()/1.e3; // mm -> m
    }
}

namespace CTPPSAlCa
{

    class RPAlignmentConstants
    {
      public:
        struct Quantities {
            Quantities() : x( 0. ), err_x( 0. ), y( 0. ), err_y( 0. ) {}
            float x, err_x, y, err_y; // in mm
        };
        RPAlignmentConstants() {
            map_[2] = Quantities();
            map_[3] = Quantities();
            map_[102] = Quantities();
            map_[103] = Quantities();
        }
        void setQuantities( unsigned int detid, const Quantities& quant ) {
            if ( map_.find( detid )==map_.end() ) return; //FIXME
            map_[detid] = quant;
        }
        const Quantities quantities( unsigned int detid ) const {
            std::map<unsigned int,Quantities>::const_iterator it = map_.find( detid );
            if ( it==map_.end() ) return Quantities();
            return it->second;
        }

      private:
        std::map<unsigned int,Quantities> map_;
    };

    struct RPCalibrationConstants {
        RPCalibrationConstants() :
            x_disp_l_n( 0. ), x_disp_l_f( 0. ),
            x_disp_r_n( 0. ), x_disp_r_f( 0. ) {}
        RPCalibrationConstants( float l_n, float l_f, float r_n, float r_f ) :
            x_disp_l_n( l_n ), x_disp_l_f( l_f ),
            x_disp_r_n( r_n ), x_disp_r_f( r_f ) {}

        // dispersions (x-axis)
        float x_disp_l_n, x_disp_l_f;
        float x_disp_r_n, x_disp_r_f;
    };

    RPCalibrationConstants getCalibrationConstants( const edm::RunNumber_t& )
    {
        RPCalibrationConstants cc;
        cc.x_disp_l_f = 9.22e-2;
        cc.x_disp_l_n = 9.26e-2;
        cc.x_disp_r_f = 5.16e-2;
        cc.x_disp_r_n = 5.81e-2;
        return cc;
    }

}

#endif
