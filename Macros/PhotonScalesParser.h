#include "TFile.h"
#include "TH2.h"
#include <iostream> //FIXME

class PhotonScalesParser
{
  public:
    PhotonScalesParser( const char* filename ) :
      file_( 0 ), hist_effmc_( 0 ) {
      file_ = new TFile( filename );
      parseFile();
    }
    ~PhotonScalesParser() {
      if ( hist_effmc_ ) delete hist_effmc_;
      if ( file_ ) delete file_;
    }

    float efficiency( float pt, float eta_sc ) const {
      if ( !hist_effmc_ ) return 0.;
      int bin_x = min( max( hist_effmc_->GetXaxis()->FindBin( eta_sc ), 1 ), hist_effmc_->GetXaxis()->GetNbins() ),
          bin_y = min( max( hist_effmc_->GetYaxis()->FindBin( pt ), 1 ), hist_effmc_->GetYaxis()->GetNbins() );
      return hist_effmc_->GetBinContent( bin_x, bin_y );
    }

  private:
    void parseFile() {
      if ( !file_ || !file_->IsOpen() ) return;
      hist_effmc_ = (TH2D*)file_->Get( "EGamma_SF2D" )->Clone();
    }

    TFile* file_;
    TH2D* hist_effmc_;
};
