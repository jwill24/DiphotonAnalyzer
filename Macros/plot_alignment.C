#include "DiphotonAnalyzer/EventAnalyzer/interface/FillNumberLUTHandler.h"
#include "DiphotonAnalyzer/EventAnalyzer/interface/AlignmentLUTHandler.h"

#include "TGraphErrors.h"
#include "TSystem.h"
#include "Canvas.h"

void plot_alignment()
{
  gSystem->Load( "libDiphotonAnalyzerEventAnalyzer.so" );

  CTPPSAlCa::FillNumberLUTHandler fill_lut( "data/fill_run_lut_v2.dat" );
  CTPPSAlCa::AlignmentLUTHandler align_lut( "data/alignment_collection_v2.out" );

  const unsigned short num_rp = 4;
  unsigned short rp_ids[num_rp] = { 2, 3, 102, 103 };
  const char* rp_names[num_rp] = { "45N", "45F", "56N", "56F" };

  vector<unsigned short> fills = fill_lut.fills();

  TGraphErrors gr_x[num_rp], gr_y[num_rp];

  unsigned short i = 0;
  for ( vector<unsigned short>::const_iterator it=fills.begin(); it!=fills.end(); it++ ) {
    CTPPSAlCa::RPAlignmentConstants ac = align_lut.getAlignmentConstants( *it );
    /*bool valid = true;
    for ( unsigned short j=0; j<num_rp; j++ ) {
      valid &= ( ac.quantities( rp_ids[j] ).x!=0. && ac.quantities( rp_ids[j] ).y!=0. );
    }
    if ( !valid ) continue;*/
    if ( ac==CTPPSAlCa::RPAlignmentConstants() ) continue;
    cout << ">> " << *it << ": " << ac << endl;
    for ( unsigned short j=0; j<num_rp; j++ ) {
      CTPPSAlCa::RPAlignmentConstants::Quantities quant = ac.quantities( rp_ids[j] );
      gr_x[j].SetPoint( i, *it, quant.x );
      gr_x[j].SetPointError( i, 0., quant.err_x );
      gr_y[j].SetPoint( i, *it, quant.y );
      gr_y[j].SetPointError( i, 0., quant.err_y );
    }
    i++;
  }

  for ( unsigned short i=0; i<num_rp; i++ ) {
    Canvas c( Form( "alignment_offsets_x_%s", rp_names[i] ), Form( "CT-PPS 2016, %s pot", rp_names[i] ) );
    gr_x[i].Draw( "ap" );
    gr_x[i].SetTitle( "LHC fill number@@x-shift (mm)" );
    gr_x[i].SetMarkerStyle( 24 );
    gr_x[i].GetYaxis()->SetRangeUser( -4.5, 0.5 );
    c.Prettify( gr_x[i].GetHistogram() );
    c.Save( "png,pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/tmp" );
  }
}
