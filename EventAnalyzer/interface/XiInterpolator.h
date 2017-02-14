#ifndef RecoCTPPS_ProtonProducer_XiInterpolator_h
#define RecoCTPPS_ProtonProducer_XiInterpolator_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/Provenance/interface/RunID.h"

#include "DiphotonAnalyzer/EventAnalyzer/interface/ProtonKinematicsUtils.h"
#include "DiphotonAnalyzer/EventAnalyzer/interface/RPAlignmentConstants.h"
#include "DiphotonAnalyzer/EventAnalyzer/interface/RPCalibrationConstants.h"

#include "TFile.h"
#include "TGraph.h"
#include "TSpline.h"

#include <iostream>

namespace ProtonUtils
{

    class XiInterpolator
    {
    public:
        XiInterpolator() :
            igLF_(0), igLN_(0), igRF_(0), igRN_(0),
            isLF_(0), isLN_(0), isRF_(0), isRN_(0)
        {}
        XiInterpolator( const char* filename ) :
            igLF_(0), igLN_(0), igRF_(0), igRN_(0),
            isLF_(0), isLN_(0), isRF_(0), isRN_(0)
        {
            loadInterpolationGraphs( filename );
        }

        ~XiInterpolator();

        void loadInterpolationGraphs( const char* filename, bool conversion=false );
        void setCalibrationConstants( const edm::RunNumber_t& run_id );
        void setAlignmentConstants( const CTPPSAlCa::RPAlignmentConstants& ac ) { align_ = ac; }

        void computeXiLinear( const TotemRPDetId& detid, const TotemRPLocalTrack& trk, float* xi_, float* err_xi_ );
        void computeXiSpline( const TotemRPDetId& detid, const TotemRPLocalTrack& trk, float* xi, float* err_xi );

    private:
        void extractSpline( const TGraph* gr_in, const char* name_out, TGraph* gr_out, TSpline3* sp_out );

        TGraph *igLF_, *igLN_, *igRF_, *igRN_;
        TSpline3 *isLF_, *isLN_, *isRF_, *isRN_;
        CTPPSAlCa::RPAlignmentConstants align_;
        CTPPSAlCa::RPCalibrationConstants calib_;
        
    };
}

#endif
