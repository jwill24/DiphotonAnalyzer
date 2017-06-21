#include "DiphotonAnalyzer/TreeProducer/interface/SelectionUtils.h"

bool
passSinglePhotonCuts( const flashgg::Photon* phot )
{
  const float abseta = fabs( phot->superCluster()->eta() );

  if ( abseta>=1.4442 and abseta<=1.566 ) return false;
  if ( phot->full5x5_r9()<=0.8 and phot->egChargedHadronIso()>=20 and phot->egChargedHadronIso()/phot->pt()>=0.3 ) return false;
  if ( phot->hadronicOverEm()>=0.08 ) return false;

  return true;
}

void
computeXiReco( const float& sqrts, const flashgg::DiPhotonCandidate* dpc, float* xi1, float* xi2 )
{
  *xi1 = ( dpc->leadingPhoton()->pt()*exp(  dpc->leadingPhoton()->eta() ) + dpc->subLeadingPhoton()->pt()*exp(  dpc->subLeadingPhoton()->eta() ) )/sqrts;
  *xi2 = ( dpc->leadingPhoton()->pt()*exp( -dpc->leadingPhoton()->eta() ) + dpc->subLeadingPhoton()->pt()*exp( -dpc->subLeadingPhoton()->eta() ) )/sqrts;
}
