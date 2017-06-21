#ifndef SelectionUtils_h
#define SelectionUtils_h

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
//#include "flashgg/DataFormats/interface/Proton.h"

bool passSinglePhotonCuts( const flashgg::Photon* phot );
void computeXiReco( const float& sqrts, const flashgg::DiPhotonCandidate* dpc, float* xi1, float* xi2 );

#endif
