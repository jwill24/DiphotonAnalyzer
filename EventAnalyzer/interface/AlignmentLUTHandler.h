#ifndef RecoCTPPS_ProtonProducer_AlignmentLUTHandler_h
#define RecoCTPPS_ProtonProducer_AlignmentLUTHandler_h

#include "DataFormats/Provenance/interface/RunID.h"
#include "DiphotonAnalyzer/EventAnalyzer/interface/RPAlignmentConstants.h"

#include <regex>
#include <fstream>

namespace CTPPSAlCa
{
    class AlignmentLUTHandler
    {
      public:
        AlignmentLUTHandler( const char* file );

        RPAlignmentConstants getAlignmentConstants( const unsigned short& fill_num ) const;
        inline bool isValid() const { return valid_; }

      private:
        // FIXME to be replaced by a DB handler!
        void loadConstants( const char* file );

        bool valid_;
        std::map<unsigned int,RPAlignmentConstants> align_map_; // fill# -> parameters
    };
}

#endif
