#include "pot_alignment.h"

void test_alignment()
{
  pot_align::load_file( "TreeProducer/data/alignment_collection_v2.out" );
  pot_align::print();

  int fill_num = 5263;

  cout << "will retrieve fill #" << fill_num << endl;
  const auto& info = pot_align::get_alignments( fill_num );
  cout << info.size() << " pots retrieved" << endl;
  for ( const auto& pot_info : info ) {
    cout << "---> pot " << pot_info.first << " -> x=" << pot_info.second.x << endl;
  }
}
