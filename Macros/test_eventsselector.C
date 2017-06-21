#include "EventsSelector.h"
#include "TSystem.h"

void
test_eventsselector()
{
  gSystem->Load( "libEventFilterUtilities.so" );
  EventsSelector ev( "/afs/cern.ch/user/l/lforthom/public/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_PPSruns_preTS2.txt" );
  cout << ev.isSelected( 273724, 1, 1 ) << std::endl; // should be out (run)
  cout << ev.isSelected( 274241, 100, 10 ) << std::endl; // should be in
  cout << ev.isSelected( 274241, 1153, 10 ) << std::endl; // should be out (ls)
}
