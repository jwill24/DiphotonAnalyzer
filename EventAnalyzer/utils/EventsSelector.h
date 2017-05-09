#ifndef EventsSelector_h
#define EventsSelector_h

#include <fstream>
#include <map>
#include <streambuf>
#include <limits.h>
#include <stdlib.h>

#include "EventFilter/Utilities/interface/reader.h"
#include "DataFormats/Provenance/interface/EventRange.h"
#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"
#include "FWCore/Utilities/interface/Exception.h"

class EventsSelector
{
  public:
    EventsSelector() {}
    EventsSelector( const char* json_file ) {
      std::ifstream file( json_file );
      std::string str( ( std::istreambuf_iterator<char>( file ) ), std::istreambuf_iterator<char>() );
      Json::Reader reader;
      Json::Value tree;

      if ( !reader.parse( str, tree ) ) throw cms::Exception( "InvalidFile" ) << "Failed to parse file " << json_file;
      if ( !tree.isObject() ) throw cms::Exception( "InvalidJSON" ) << "Invalid format!";

      Json::Value::Members runs = tree.getMemberNames(); //vector<string>
      for ( Json::Value::Members::const_iterator it_runs=runs.begin(); it_runs!=runs.end(); it_runs++ ) {
	const edm::RunNumber_t run = atoi( it_runs->c_str() );
	Json::Value lumi_ranges = tree[*it_runs];
	for ( Json::Value::iterator rng=lumi_ranges.begin(); rng!=lumi_ranges.end(); rng++ ) {
	  if ( ( *rng ).size()!=2 ) throw cms::Exception( "InvalidJSON" ) << "Invalid format!";
	  const edm::LuminosityBlockNumber_t ls_beg = ( *rng )[0u].asUInt(),
                                             ls_end = ( *rng )[1u].asUInt();
	  ranges_.insert( std::make_pair( run, edm::EventRange( run, ls_beg, 1, run, ls_end, UINT_MAX ) ) );
	}
      }
    }
    ~EventsSelector() {}

    unsigned short size() const { return ranges_.size(); }
    bool isSelected( unsigned int run, unsigned int ls, unsigned long long evt ) const {
      if ( ranges_.count( run )==0 ) return false;
      const edm::EventID evt_id( run, ls, evt );
      for ( map_t::const_iterator it=ranges_.lower_bound( run ); it!=ranges_.upper_bound( run ); it++ ) {
	if ( evt_id>=it->second.startEventID() && evt_id<=it->second.endEventID() ) return true;
      }
      return false;
    }

  private:
    typedef std::map<edm::RunNumber_t,edm::EventRange> map_t;
    map_t ranges_;
};

#endif
