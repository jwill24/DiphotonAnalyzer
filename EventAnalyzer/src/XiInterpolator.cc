#include "DiphotonAnalyzer/EventAnalyzer/interface/XiInterpolator.h"

namespace ProtonUtils
{
    void
    XiInterpolator::loadInterpolationGraphs( const char* filename, bool conversion )
    {
        TFile f( filename );
        if ( !f.IsOpen() ) {
            edm::LogError("XiInterpolator") << "Failed to load the interpolation graphs file";
            return;
        }

        if ( conversion ) {
            // in case one works with Frici's input
            extractSpline( (TGraph*)f.Get("XRPH_C6R5_B1"), "x_to_xi_R_1_N", igRN_, isRN_ );
            extractSpline( (TGraph*)f.Get("XRPH_D6R5_B1"), "x_to_xi_R_1_F", igRF_, isRF_ );
            extractSpline( (TGraph*)f.Get("XRPH_C6L5_B2"), "x_to_xi_L_1_N", igLN_, isLN_ );
            extractSpline( (TGraph*)f.Get("XRPH_D6L5_B2"), "x_to_xi_L_1_F", igLF_, isLF_ );
            return;
        }
        // already "processed" curves
        igRN_ = (TGraph*)f.Get("g_x_to_xi_R_1_N");
        igRF_ = (TGraph*)f.Get("g_x_to_xi_R_1_F");
        igLN_ = (TGraph*)f.Get("g_x_to_xi_L_1_N");
        igLF_ = (TGraph*)f.Get("g_x_to_xi_L_1_F");
        isRN_ = (TSpline3*)f.Get("s_x_to_xi_R_1_N");
        isRF_ = (TSpline3*)f.Get("s_x_to_xi_R_1_F");
        isLN_ = (TSpline3*)f.Get("s_x_to_xi_L_1_N");
        isLF_ = (TSpline3*)f.Get("s_x_to_xi_L_1_F");

        edm::LogInfo("XiInterpolator") << "Interpolation graphs successfully loaded from file " << filename;
    }

    void
    XiInterpolator::setCalibrationConstants( const edm::RunNumber_t& run_id )
    {
        calib_ = CTPPSAlCa::getCalibrationConstants( run_id );
        edm::LogInfo("ProtonReco")
            << "Calibration constants loaded: " << calib_.x_disp_l_f << ":" << calib_.x_disp_l_n << " / "
                                                << calib_.x_disp_r_f << ":" << calib_.x_disp_r_n;
    }

    void
    XiInterpolator::computeXiLinear( const TotemRPDetId& detid, const TotemRPLocalTrack& trk, float* xi_, float* err_xi_ )
    {
        *xi_ = *err_xi_ = 0.;

        if ( !trk.isValid() ) return;

        float dx_n, dx_f;
        const CTPPSAlCa::RPAlignmentConstants::Quantities ac = align_.quantities( detid.rawId() );
        switch ( detid.arm() ) { // 0 = left, 1 = right
            case 0: {
                dx_n = calib_.x_disp_l_n;
                dx_f = calib_.x_disp_l_f;
            } break;
            case 1: {
                dx_n = calib_.x_disp_r_n;
                dx_f = calib_.x_disp_r_f;
            } break;
            default: return;
        }

        const float de_x = 0.2e-3/*m*/, de_rel_dx = 0.1;
        const float x_corr = ( trk.getX0() + ac.x ) * 1.e-3;

        if ( detid.romanPot()==3 ) { // far pot //FIXME
            *xi_ = x_corr / dx_f;
            *err_xi_ = std::sqrt( std::pow( de_x/dx_f, 2 )
                                + std::pow( de_rel_dx * (*xi_), 2 ) );
            return;
        }
        if ( detid.romanPot()==2 ) { // near pot //FIXME
            *xi_ = x_corr / dx_n;
            *err_xi_ = std::sqrt( std::pow( de_x/dx_n, 2 )
                                + std::pow( de_rel_dx * (*xi_), 2 ) );
            return;
        }
    }

    void
    XiInterpolator::computeXiSpline( const TotemRPDetId& detid, const TotemRPLocalTrack& trk, float* xi, float* err_xi )
    {
        *xi = *err_xi = 0.;

        if ( !trk.isValid() ) return;

        TSpline3 *interp_near = 0, *interp_far = 0;
        const CTPPSAlCa::RPAlignmentConstants::Quantities ac = align_.quantities( detid.rawId() );
        switch ( detid.arm() ) { // 0 = left, 1 = right
            case 0: {
                interp_far = isLF_;
                interp_near = isLN_;
            } break;
            case 1: {
                interp_far = isRF_;
                interp_near = isRN_;
            } break;
            default: return;
        }
        if ( ( !interp_far ) or ( !interp_near ) ) return;

        const float de_x = 0.4e-3/*m*/, de_rel_dx = 0.1;
        const float x_corr = ( trk.getX0() + ac.x ) * 1.e-3;

        if ( detid.romanPot()==3 ) { // far pot //FIXME
            *xi = interp_far->Eval( x_corr );
            const float de_xi = interp_far->Eval( x_corr + de_x ) - ( *xi );
            *err_xi = std::sqrt( std::pow( de_xi, 2 )
                               + std::pow( de_rel_dx * (*xi), 2 ) );
            return;
        }
        if ( detid.romanPot()==2 ) { // near pot //FIXME
            *xi = interp_near->Eval( x_corr );
            const float de_xi = interp_near->Eval( x_corr + de_x ) - ( *xi );
            *err_xi = std::sqrt( std::pow( de_xi, 2 )
                               + std::pow( de_rel_dx * (*xi), 2 ) );
            return;
        }
    }

    void
    XiInterpolator::extractSpline( const TGraph* gr_in, const char* name_out, TGraph* gr_out, TSpline3* sp_out )
    {
        if ( !gr_in ) return;

        const double offset = gr_in->GetX()[0];
        if ( gr_out ) delete gr_out;
        gr_out = new TGraph( Form( "g_%s", name_out ), gr_in->GetTitle() );

        double x, y;
        for ( int i=0; i<gr_in->GetN(); i++ ) {
            gr_in->GetPoint( i, x, y );
            gr_out->SetPoint( i, x-offset, -y );
        }
        if ( sp_out ) delete sp_out;
        sp_out = new TSpline3( "", gr_out->GetX(), gr_out->GetY(), gr_out->GetN() );
        sp_out->SetName( Form( "s_%s", name_out ) );
    }
}
