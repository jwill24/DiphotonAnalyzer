#ifndef DiphotonAnalyzer_Macros_xi_reconstruction_h
#define DiphotonAnalyzer_Macros_xi_reconstruction_h

#include "TSpline.h"

namespace xi_reco
{
  TSpline* sp[4];

  void
  load_file( const char* filename )
  {
    TFile* file = new TFile( filename );
    sp[0] = dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_L_1_N" )->Clone() ); //45N
    sp[1] = dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_L_1_F" )->Clone() ); //45F
    sp[2] = dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_R_1_N" )->Clone() ); //56N
    sp[3] = dynamic_cast<TSpline*>( file->Get( "s_x_to_xi_R_1_F" )->Clone() ); //56F
    delete file;
  }

  TSpline*
  get_spline( unsigned int arm, unsigned int pot )
  {
    if ( arm==0 && pot==2 ) return sp[0];
    if ( arm==0 && pot==3 ) return sp[1];
    if ( arm==1 && pot==2 ) return sp[2];
    if ( arm==1 && pot==3 ) return sp[3];
    return 0;
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
}

#endif
