#ifndef RecoCTPPS_ProtonProducer_RPAlignmentConstants_h
#define RecoCTPPS_ProtonProducer_RPAlignmentConstants_h

#include <map>

namespace CTPPSAlCa
{
    class RPAlignmentConstants
    {
      public:
        typedef struct Quantities {
            Quantities() : x( 0. ), err_x( 0. ), y( 0. ), err_y( 0. ) {}
            float x, err_x, y, err_y; // in mm
        } Quantities;
        RPAlignmentConstants();
        void setQuantities( unsigned int detid, const Quantities& quant );
        const Quantities quantities( unsigned int detid ) const;

      private:
        std::map<unsigned int,Quantities> map_;
    };
}

#endif
