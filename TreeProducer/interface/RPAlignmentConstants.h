#ifndef RecoCTPPS_ProtonProducer_RPAlignmentConstants_h
#define RecoCTPPS_ProtonProducer_RPAlignmentConstants_h

#include <iostream>
#include <map>

namespace CTPPSAlCa
{
  class RPAlignmentConstants
  {
    public:
      class Quantities
      {
        public:
          Quantities() : x( 0. ), err_x( 0. ), y( 0. ), err_y( 0. ) {}
          Quantities( float x, float y, float err_x=0., float err_y=0. ) : x( x ), err_x( err_x ), y( y ), err_y( err_y ) {}

          bool operator==( const Quantities& rhs ) const {
            return ( x==rhs.x && err_x==rhs.err_x
                  && y==rhs.y && err_y==rhs.err_y );
          }
          void operator=( const Quantities& rhs ) {
            x = rhs.x; err_x = rhs.err_x;
            y = rhs.y; err_y = rhs.err_y;
          }
          friend std::ostream& operator<<( std::ostream&, const Quantities& );

        public:
          float x, err_x, y, err_y;
      };
      typedef std::map<unsigned short, Quantities> Map;

    public:
      RPAlignmentConstants();
      friend std::ostream& operator<<( std::ostream&, const RPAlignmentConstants& );
      bool operator==( const RPAlignmentConstants& ) const;
      void operator=( const RPAlignmentConstants& );

      Map::const_iterator begin() const { return map_.begin(); }
      Map::const_iterator end() const { return map_.end(); }

      void setQuantities( unsigned short detid, const Quantities& quant );
      const Quantities quantities( unsigned short detid ) const;

    private:
      Map map_;
  };
}

#endif
