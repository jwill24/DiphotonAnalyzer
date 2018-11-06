#ifndef DiphotonAnalyzer_Macros_xi_reconstruction_2017_h
#define DiphotonAnalyzer_Macros_xi_reconstruction_2017_h

#include "TSpline.h"
#include "TFile.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <map>

namespace xi_reco
{
  struct dispersion_t
  {
  dispersion_t( double value = 0., double error = 0. ) :
    value( value ), err_value( error ) {}
    double value, err_value;
  };
  typedef std::map<unsigned short,dispersion_t> dispersions_t;
  class xi_reconstruction_2017
  {
  public:
    xi_reconstruction_2017() {}
    ~xi_reconstruction_2017() {}


    /* Naming convention:
       B1 = SECTOR 56 (-z in the CMS coordinate system)
       B2 = SECTOR 45 (+z in the CMS coordinate system)
    */


    //Open dispersion files    
    TFile *f_inb1_120 = TFile::Open("TreeProducer/data/xi_as_a_function_of_x_graph_b1_120_murad.root");
    TFile *f_inb1_130 = TFile::Open("TreeProducer/data/xi_as_a_function_of_x_graph_b1_130_murad.root");
    TFile *f_inb1_140 = TFile::Open("TreeProducer/data/xi_as_a_function_of_x_graph_b1_140_murad.root");
    TFile *f_inb2_120 = TFile::Open("TreeProducer/data/xi_as_a_function_of_x_graph_b2_120_murad.root");
    TFile *f_inb2_130 = TFile::Open("TreeProducer/data/xi_as_a_function_of_x_graph_b2_130_murad.root");
    TFile *f_inb2_140 = TFile::Open("TreeProducer/data/xi_as_a_function_of_x_graph_b2_140_murad.root");
    
    TGraph *xtoxib1_120 = (TGraph *) f_inb1_120->Get("XRPH_D6R5_B1");
    TGraph *xtoxib1_130 = (TGraph *) f_inb1_130->Get("XRPH_D6L5_B1");
    TGraph *xtoxib1_140 = (TGraph *) f_inb1_140->Get("XRPH_D6R5_B1");
    TGraph *xtoxib2_120 = (TGraph *) f_inb2_120->Get("XRPH_D6L5_B2");
    TGraph *xtoxib2_130 = (TGraph *) f_inb2_130->Get("XRPH_D6L5_B2");
    TGraph *xtoxib2_140 = (TGraph *) f_inb2_140->Get("XRPH_D6L5_B2");
    

    //Load my complementary dispersion file
    void feedDispersions( const char* filename ) {
      std::ifstream file( filename );
      if ( !file.is_open() )
	throw std::runtime_error( std::string( "Failed to open ")+filename );
      std::string line;
      while ( !file.eof() ) {
	std::getline( file, line );
	if ( line[0] == '#' ) continue; // strip comments
	std::istringstream ss( line );
	std::string s_xangle, s_pot, s_disp, s_err_disp;
	while ( ss >> s_xangle >> s_pot >> s_disp >> s_err_disp ) {
	  const float xangle = std::stod( s_xangle );
	  const unsigned short pot = std::stoi( s_pot );
	  const double disp = std::stod( s_disp ), err_disp = std::stod( s_err_disp );
	  disp_vs_xa_[xangle][pot] = dispersion_t( disp, err_disp );
	}
      }
    }
    
    void reconstruct( double xangle, unsigned short pot, double x_aligned, double& xi, double& xi_err ) const {
      if ( pot == 103 || pot == 123 ) {
	/*
	if ( xangle == 120 ) xi = -0.01* xtoxib1_120->Eval( x_aligned );
	else if ( xangle == 130 ) xi = -0.01* xtoxib1_130->Eval( x_aligned );
	else if ( xangle == 140 ) xi = -0.01* xtoxib1_140->Eval( x_aligned );
	else if ( xangle < 120 ) xi = fabs( xtoxib1_120->Eval( x_aligned ) + ( ( ( xangle - 120. )/( 120. - 140. ) )*( xtoxib1_140->Eval( x_aligned )-xtoxib1_120->Eval( x_aligned ) ) ) );
	else if ( xangle > 120 && xangle < 150 ) xi = fabs( xtoxib1_120->Eval( x_aligned ) + ( ( ( 120. - xangle )/( 120. - 140. ) )*( xtoxib1_140->Eval( x_aligned )-xtoxib1_120->Eval( x_aligned ) ) ) );
	*/
	//else {
	const auto& disp_vs_pot = disp_vs_xa_.at( xangle );
	const dispersion_t disp = disp_vs_pot.at( pot );
	xi = x_aligned/( -disp.value );
	  //}
      }
      
      if ( pot == 3 || pot == 23 ) {
	/*
	if ( xangle == 120 ) xi = -0.01* xtoxib2_120->Eval( x_aligned );
	else if ( xangle == 130 ) xi = -0.01* xtoxib2_130->Eval( x_aligned );
	else if ( xangle == 140 ) xi = -0.01* xtoxib2_140->Eval( x_aligned );
	else if ( xangle < 120 ) xi = fabs( xtoxib1_120->Eval( x_aligned ) + ( ( ( xangle - 120. )/( 120. - 140. ) )*( xtoxib1_140->Eval( x_aligned )-xtoxib1_120->Eval( x_aligned ) ) ) );
	else if ( xangle > 120 && xangle < 150 ) xi = fabs( xtoxib2_120->Eval( x_aligned ) + ( ( ( 120. - xangle )/( 120. - 140. ) )*( xtoxib2_140->Eval( x_aligned )-xtoxib2_120->Eval( x_aligned ) ) ) );
	*/
	//else {
	const auto& disp_vs_pot = disp_vs_xa_.at( xangle );
	const dispersion_t disp = disp_vs_pot.at( pot );
	xi = x_aligned/( -disp.value );
	//}
      }	  
    }
    
    
 private:
    std::map<float,dispersions_t> disp_vs_xa_;
  };
}

#endif


    
