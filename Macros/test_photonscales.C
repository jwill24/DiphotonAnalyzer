#include "PhotonScalesParser.h"
#include <iostream>

void test_photonscales()
{
  PhotonScalesParser parser( "egammaEffi.txt_EGM2D.root" );
  cout << ">>> " << parser.efficiency( 550., 2.6 ) << endl;
}
