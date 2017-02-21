#include "DiphotonAnalyzer/EventAnalyzer/interface/RPAlignmentConstants.h"

namespace CTPPSAlCa
{
    RPAlignmentConstants::RPAlignmentConstants()
    {
        map_[2] = Quantities();
        map_[3] = Quantities();
        map_[102] = Quantities();
        map_[103] = Quantities();
    }

    void
    RPAlignmentConstants::setQuantities( unsigned int detid, const Quantities& quant )
    {
        if ( map_.find( detid )==map_.end() ) return; //FIXME
        map_[detid] = quant;
    }

    const RPAlignmentConstants::Quantities
    RPAlignmentConstants::quantities( unsigned int detid ) const
    {
        std::map<unsigned int,Quantities>::const_iterator it = map_.find( detid );
        if ( it==map_.end() ) return Quantities();
        return it->second;
    }
}
