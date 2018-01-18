#ifndef DiphotonAnalyzer_Macros_xi_reconstruction_h
#define DiphotonAnalyzer_Macros_xi_reconstruction_h

#include "TSpline.h"
#include <map>

namespace xi_reco
{
  std::map<unsigned short,TSpline*> splines;
  std::map<unsigned short,std::pair<double,double> > dispersions; // RP -> ( disp, err_disp )

  void set_dispersions( const std::map<unsigned short,std::pair<double,double> >& disp ) { dispersions = disp; }

  void
  load_file( const char* filename )
  {
    auto file = TFile::Open( filename );
    splines = std::map<unsigned short,TSpline*>( {
      { 2, dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_L_1_N" )->Clone() ) }, //45N
      { 3, dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_L_1_F" )->Clone() ) }, //45F
      { 102, dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_R_1_N" )->Clone() ) }, //56N
      { 103, dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_R_1_F" )->Clone() ) } //56F
    } );
    delete file;
  }

  TSpline*
  get_spline( unsigned int arm, unsigned int pot )
  {
    auto it_sp = splines.find( 100*arm+pot );
    if ( it_sp == splines.end() ) return 0;
    return it_sp->second;
  }

  void
  reconstruct( double x, unsigned int arm, unsigned int pot, double& xi, double& xi_err )
  {
    xi = xi_err = 0.;
    TSpline* sp = get_spline( arm, pot );
    if ( !sp ) return;
    xi = sp->Eval( x );
    // old
    /*double de_xi = sp->Eval( x+0.4e-3 )-xi, de_dx = 0.1*xi;
    xi_err = sqrt( de_xi*de_xi+de_dx*de_dx );*/
    // determine uncertainty of xi
    const double si_x_alignment = 150.0e-6; // in m, alignment uncertainty
    const double si_x_neglected_angle = 150.0e-6; // in m, to (approximately) account for the neglected angular term in proton transport
    const double si_rel_D = 0.055; // 1, relative uncertainty of dispersion
    const double si_x = sqrt( si_x_alignment*si_x_alignment + si_x_neglected_angle*si_x_neglected_angle );
    const double si_xi_from_x = sp->Eval( x+si_x )-xi;
    const double si_xi_from_D_x = si_rel_D * xi;
    xi_err = sqrt( si_xi_from_x*si_xi_from_x + si_xi_from_D_x*si_xi_from_D_x );
  }

  void
  reconstruct_lin( double x, unsigned int arm, unsigned int pot, double& xi, double& xi_err )
  {
    xi = xi_err = 0;
    if ( dispersions.count( 100*arm+pot ) == 0 ) return;

    const auto disp = dispersions.at( 100*arm+pot ); // ( disp, err_disp )
    const double de_x = 150.e-6; // alignment uncertainty

    xi = x/disp.first;
    err_xi = sqrt( pow( de_x/disp.first, 2 )+pow( disp.second*xi, 2 ) );
  }
}

#endif
