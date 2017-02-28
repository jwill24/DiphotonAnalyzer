#include "DiphotonAnalyzer/EventAnalyzer/interface/XiInterpolator.h"

namespace ProtonUtils
{
  XiInterpolator::XiInterpolator() :
    igLF_(0), igLN_(0), igRF_(0), igRN_(0),
    isLF_(0), isLN_(0), isRF_(0), isRN_(0)
  {}

  XiInterpolator::XiInterpolator( const char* filename ) :
    igLF_(0), igLN_(0), igRF_(0), igRN_(0),
    isLF_(0), isLN_(0), isRF_(0), isRN_(0)
  {
    loadInterpolationGraphs( filename );
  }

  XiInterpolator::~XiInterpolator()
  {}

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
  XiInterpolator::computeXiLinear( const TotemRPDetId& detid, const TotemRPLocalTrack& trk, float& xi, float& err_xi )
  {
    xi = err_xi = 0.;

    if ( !trk.isValid() ) return;

    // retrieve the alignment parameters
    const CTPPSAlCa::RPAlignmentConstants::Quantities ac = align_.quantities( detid.rawId() );

    // retrieve the proper dispersion constants
    float dx_n, dx_f;
    switch ( detid.arm() ) { // 0 = sector 45, 1 = sector 56
      case 0: { dx_n = calib_.x_disp_l_n; dx_f = calib_.x_disp_l_f; } break;
      case 1: { dx_n = calib_.x_disp_r_n; dx_f = calib_.x_disp_r_f; } break;
      default: return;
    }

    const float de_x = 0.2e-3, // m
                de_rel_dx = 0.1;

    // apply the alignment
    const float x_corr = ( trk.getX0() + ac.x ) * 1.e-3;

    if ( detid.romanPot()==3 ) { // far pot
      xi = x_corr / dx_f;
      err_xi = std::sqrt( std::pow( de_x/dx_f, 2 )
                        + std::pow( de_rel_dx * xi, 2 ) );
    }
    if ( detid.romanPot()==2 ) { // near pot
      xi = x_corr / dx_n;
      err_xi = std::sqrt( std::pow( de_x/dx_n, 2 )
                        + std::pow( de_rel_dx * xi, 2 ) );
    }
  }

  void
  XiInterpolator::computeXiSpline( const TotemRPDetId& detid, const TotemRPLocalTrack& trk, float& xi, float& err_xi )
  {
    xi = err_xi = 0.;

    if ( !trk.isValid() ) return;

    // retrieve the alignment parameters
    const CTPPSAlCa::RPAlignmentConstants::Quantities ac = align_.quantities( detid.rawId() );

    // retrieve the proper interpolation curve
    TSpline3 *interp = 0;
    switch ( detid.arm() ) { // 0 = sector 45, 1 = sector 56
      case 0: {
        if ( detid.romanPot()==3 ) interp = isLF_;
        if ( detid.romanPot()==2 ) interp = isLN_;
      } break;
      case 1: {
        if ( detid.romanPot()==3 ) interp = isRF_;
        if ( detid.romanPot()==2 ) interp = isRN_;
      } break;
      default: return;
    }
    if ( !interp ) return;

    const float de_x = 0.4e-3, // m
                de_rel_dx = 0.1;

    // apply the alignment
    const float x_corr = ( trk.getX0() + ac.x ) * 1.e-3;

    xi = interp->Eval( x_corr );

    const float de_xi = interp->Eval( x_corr + de_x ) - xi;
    err_xi = std::sqrt( std::pow( de_xi, 2 )
                      + std::pow( de_rel_dx * xi, 2 ) );
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
