#include "DiphotonAnalyzer/EventAnalyzer/interface/AlignmentLUTHandler.h"

namespace CTPPSAlCa
{
    AlignmentLUTHandler::AlignmentLUTHandler( const char* file ) :
        valid_( false )
    {
        // load the alignment constants from external file
        loadConstants( file );
        valid_ = true;
    }

    RPAlignmentConstants
    AlignmentLUTHandler::getAlignmentConstants( const unsigned short& fill_num ) const
    {
        if ( !valid_ or !align_map_.count( fill_num ) ) return RPAlignmentConstants();
        return align_map_.at( fill_num );
    }

    void
    AlignmentLUTHandler::loadConstants( const char* file )
    {
        std::ifstream f( file );
        if ( !f.is_open() ) return;
        std::string ss;
        //bool has_margin = false, x_method = false;
        unsigned short fill_num = 0;

        //std::regex rgx( ".*id=(\\d+),sh_x=(-?[0-9]\\d*(\\.\\d+)?).*" );
        std::regex rgx( ".*id=(\\d+),sh_x=([0-9-\\.]+)(?:,sh_x_unc=)?([0-9-\\.]+)?(?:,sh_y=)?([0-9-\\.]+)?(?:,sh_y_unc=)?([0-9-\\.]+)?" );
        std::smatch match;

        RPAlignmentConstants ac;

        while ( f >> ss ) {
            size_t pos1 = ss.find( "[" );
            if ( pos1!=std::string::npos ) { // new fill info
                align_map_.insert( std::pair<unsigned int,RPAlignmentConstants>( fill_num, ac ) );
                ac = RPAlignmentConstants(); // reset the constants list
                /*has_margin = ( ss.find( "no_margin" )==std::string::npos ); // .5mm extra margin
                x_method = ( ss.find( "method x" )!=std::string::npos ); // computation method*/
                //size_t pos2 = ss.find( "fill_" ), pos3 = ss.find( "/method" );
                size_t pos2 = ss.find( "fill_" ), pos3 = ss.find( "/201" ); //FIXME use a regex here too!
                fill_num = atoi( ss.substr( pos2+5, pos3 ).c_str() );
                //std::cout << "-----> fill_num=" << fill_num << ", margin:" << has_margin << ",x_method:" << x_method << std::endl;
                continue;
            }
            const std::string s = ss;
            if ( std::regex_search( s.begin(), s.end(), match, rgx ) ) {
                const unsigned int rp_id = atoi( match[1].str().c_str() );
                const float align_x = atof( match[2].str().c_str() );
                float err_align_x = 0., align_y = 0., err_align_y = 0.;
                if ( match.size()>3 ) {
                    err_align_x = atof( match[3].str().c_str() );
                    align_y = atof( match[4].str().c_str() );
                    err_align_y = atof( match[5].str().c_str() );
                }
                //const unsigned short arm = rp_id/100, pot = rp_id%100;
                RPAlignmentConstants::Quantities quant;
                quant.x = align_x;
                quant.err_x = err_align_x;
                quant.y = align_y;
                quant.err_y = err_align_y;
                ac.setQuantities( rp_id, quant );
            }
        }
    }
}
