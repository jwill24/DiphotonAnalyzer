#include "DiphotonAnalyzer/TreeProducer/interface/RPAlignmentConstants.h"

namespace CTPPSAlCa
{
  RPAlignmentConstants::RPAlignmentConstants()
  {
    map_[2] = Quantities();
    map_[3] = Quantities();
    map_[102] = Quantities();
    map_[103] = Quantities();
  }

  bool
  RPAlignmentConstants::operator==( const RPAlignmentConstants& rhs ) const
  {
    return ( map_==rhs.map_ );
  }

  void
  RPAlignmentConstants::operator=( const RPAlignmentConstants& rhs )
  {
    map_ = rhs.map_;
  }

  void
  RPAlignmentConstants::setQuantities( unsigned short detid, const Quantities& quant )
  {
    if ( map_.find( detid )==map_.end() ) return; //FIXME
    map_[detid] = quant;
  }

  const RPAlignmentConstants::Quantities
  RPAlignmentConstants::quantities( unsigned short detid ) const
  {
    Map::const_iterator it = map_.find( detid );
    if ( it!=map_.end() ) return it->second;
    std::cerr << "WARNING failed to retrieve quantities for RP with DetId " << detid << std::endl;
    return Quantities();
  }

  std::ostream&
  operator<<( std::ostream& os, const RPAlignmentConstants& align )
  {
    for ( RPAlignmentConstants::Map::const_iterator it=align.map_.begin(); it!=align.map_.end(); ++it ) {
      os << "DetId = " << it->first << ": " << it->second << std::endl;
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const RPAlignmentConstants::Quantities& quant )
  {
    return os << "x-shift: " << quant.x << " +/- " << quant.err_x << "\t"
              << "y-shift: " << quant.y << " +/- " << quant.err_y;
  }
}
