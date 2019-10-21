/**
  The original version of this file is taken from the MiniZinc source code
  www.minizinc.org
  */

/*
 *  Main authors:
 *     Gleb Belov <gleb.belov@monash.edu>
 */

#include <map>
#include <sstream>
#include <cmath>

#include "mip_seccutgen.h"
#include "min_cut.h"
#include "utils.hh"

std::string MIP_SECCutGen::validate() const {
  std::ostringstream oss;
  /// Check that diagonal flows are 0
  for (int i=0; i<nN; ++i)
    if (pMIP->colUB[varXij[i*nN+i]] > 0.0)
      oss << "SECutGen with " << nN
          << " cities: diagonal flow " << (i+1)
          << " has UB=" << pMIP->colUB[varXij[i*nN+i]]
          << "\n";
  return oss.str();
}

void MIP_SECCutGen::generate(const MIP_wrapper::Output& slvOut, MIP_wrapper::CutInput& cutsIn) {
  assert( pMIP );
  /// Extract graph, converting to undirected
  typedef std::map< std::pair< int, int >, double > TMapFlow;
  TMapFlow mapFlow;
  for ( int i=0; i<nN; ++i ) {
    for ( int j=0; j<nN; ++j ) {
      const double xij = slvOut.x[ varXij[ nN*i + j ] ];
      if ( i==j )
        MZN_ASSERT_HARD_MSG( 1e-4 > std::fabs(xij), "circuit: X[" << (i+1) << ", " << (j+1) << "]==" << xij );
      MZN_ASSERT_HARD_MSG( -1e-4 < xij && 1.0+1e-4 > xij,           // adjusted from 1e-6 to 1e-4 for CBC. 7.8.19
                           "circuit: X[" << (i+1) << ", " << (j+1) << "]==" << xij );
      if ( 1e-4 <= xij ) {
        mapFlow[ std::make_pair( std::min(i,j), std::max(i,j) ) ] += xij;
      }
    }
  }
  /// Invoking Min Cut
//   cerr << "  MIN CUT... " << flush;
  Algorithms::MinCut mc;
  mc.nNodes = nN;
  mc.edges.reserve( mapFlow.size() );
  mc.weights.reserve( mapFlow.size() );
  for ( const auto& mf: mapFlow ) {
    mc.edges.push_back( mf.first );
    mc.weights.push_back( mf.second );
  }
  mc.solve();
  /// Check if violation
  if ( mc.wMinCut <= 1.999 ) {
    MIP_wrapper::CutDef cut( MIP_wrapper::GQ, MIP_wrapper::MaskConsType_Lazy | MIP_wrapper::MaskConsType_Usercut );
    cut.rhs = 1.0;
    int nCutSize=0;
    constexpr int nElemPrint = 20;
    // cerr << "  CUT: [ ";
    for ( int i=0; i<nN; ++i )
    if ( mc.parities[i] ) {
      ++nCutSize;
      //if ( nCutSize<=nElemPrint )
      //  cerr << (i+1) << ", ";
      //else if ( nCutSize==nElemPrint+1 )
      //  cerr << "...";
      for ( int j=0; j<nN; ++j )
      if ( !mc.parities[j] ) {
        cut.addVar( varXij[ nN*i + j ], 1.0 );
      }
    }
    // cerr << "]. " << flush;
    double dViol = cut.computeViol( slvOut.x, slvOut.nCols );
    if ( dViol > 0.0001 ) {   // ?? PARAM?  TODO. See also min cut value required
      cutsIn.push_back( cut );
      std::cerr << "  SEC cut: violation=" << dViol
        << "  N NODES: " << nN
        << "  |X|: : " << nCutSize
        << std::endl;
    } else {
      MZN_ASSERT_HARD_MSG( 0, "  SEC cut: N nodes = " << nN << ": violation = " << dViol
        << ": too small compared to the min-cut value " << (2.0-mc.wMinCut)
      );
    }
  }
}

void MIP_SECCutGen::print(std::ostream&) {
}
