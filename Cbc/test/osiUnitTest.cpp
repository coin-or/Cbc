// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"

#include "OsiConfig.h"

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>

#include "OsiUnitTests.hpp"

#include "CoinError.hpp"

#include "OsiCbcSolverInterface.hpp"

namespace {

// Display message on stdout and stderr. Flush cout buffer before printing the
// message, so that output comes out in order in spite of buffered cout.

void testingMessage( const char * const msg )
{
  std::cout.flush() ;
  std::cerr <<msg;
  //cout <<endl <<"*****************************************"
  //     <<endl <<msg <<endl;
}


/*
  Utility routine to process command line parameters. An unrecognised parameter
  will trigger the help message and a return value of false.
  
  This should be replaced with the one of the standard CoinUtils parameter
  mechanisms.
*/
bool processParameters (int argc, const char **argv,
			std::map<std::string,std::string> &parms)

{
  assert(argc >= 1);
  assert(argv != NULL);
/*
  Initialise the parameter keywords.
*/
  std::set<std::string> definedKeyWords;
  definedKeyWords.insert("-cerr2cout");
  definedKeyWords.insert("-mpsDir");
  definedKeyWords.insert("-netlibDir");
  definedKeyWords.insert("-testOsiSolverInterface");
  definedKeyWords.insert("-nobuf");
/*
  Set default values for data directories.
*/
  const char dirsep =  CoinFindDirSeparator() ;
  std::string pathTmp ;

  pathTmp = ".." ;
  pathTmp += dirsep ;
  pathTmp += ".." ;
  pathTmp += dirsep ;
  pathTmp += "Data" ;
  pathTmp += dirsep ;

  parms["-mpsDir"] = pathTmp + "Sample"  ;
  parms["-netlibDir"] = pathTmp + "Netlib" ;

/*
  Read the command line parameters and fill a map of parameter keys and
  associated data. The parser allows for parameters which are only a keyword,
  or parameters of the form keyword=value (no spaces).
*/
  for (int i = 1 ; i < argc ; i++)
  { std::string parm(argv[i]) ;
    std::string key,value ;
    std::string::size_type eqPos = parm.find('=');

    if (eqPos == std::string::npos)
    { key = parm ; }
    else
    { key = parm.substr(0,eqPos) ;
      value = parm.substr(eqPos+1) ; }
/*
  Is the specifed key valid?
*/
    if (definedKeyWords.find(key) == definedKeyWords.end())
    { std::cerr << "Undefined parameter \"" << key << "\"." << std::endl ;
      std::cerr
	<< "Usage: " << argv[0]
	<< " [-nobuf] [-mpsDir=V1] [-netlibDir=V2] "
        << "[-testOsiSolverInterface] " << std::endl ;
      std::cerr << "  where:" << std::endl ;
      std::cerr
	<< "    "
	<< "-cerr2cout: redirect cerr to cout; sometimes useful." << std::endl
	<< "\t" << "to synchronise cout & cerr." << std::endl ;
      std::cerr
	<< "    "
	<< "-mpsDir: directory containing mps test files." << std::endl
        << "\t" << "Default value V1=\"../../Data/Sample\"" << std::endl ;
      std::cerr
	<< "    "
	<< "-netlibDir: directory containing netlib files." << std::endl
        << "\t" << "Default value V2=\"../../Data/Netlib\"" << std::endl ;
      std::cerr
	<< "    "
	<< "-testOsiSolverInterface: "
        << "run each OSI on the netlib problem set." << std::endl
	<< "\t"
	<< "Default is to not run the netlib problem set." << std::endl ;
      std::cerr
	<< "    "
        << "-nobuf: use unbuffered output." << std::endl
	<< "\t" << "Default is buffered output." << std::endl ;
      
      return (false) ; }
/*
  Valid keyword; stash the value for later reference.
*/
    parms[key]=value ; }
/*
  Tack the directory separator onto the data directories so we don't have to
  worry about it later.
*/
  parms["-mpsDir"] += dirsep ;
  parms["-netlibDir"] += dirsep ;
/*
  Did the user request unbuffered i/o? It seems we need to go after this
  through stdio --- using pubsetbuf(0,0) on the C++ streams has no
  discernible affect. Nor, for that matter, did setting the unitbuf flag on
  the streams. Why? At a guess, sync_with_stdio connects the streams to the
  stdio buffers, and the C++ side isn't programmed to change them?
*/
  if (parms.find("-nobuf") != parms.end())
  { // std::streambuf *coutBuf, *cerrBuf ;
    // coutBuf = std::cout.rdbuf() ;
    // coutBuf->pubsetbuf(0,0) ;
    // cerrBuf = std::cerr.rdbuf() ;
    // cerrBuf->pubsetbuf(0,0) ;
    setbuf(stderr,0) ;
    setbuf(stdout,0) ; }
/*
  Did the user request a redirect for cerr? This must occur before any i/o is
  performed.
*/
  if (parms.find("-cerr2cout") != parms.end())
  { std::cerr.rdbuf(std::cout.rdbuf()) ; }

  return (true) ; }


}	// end file-local namespace



//----------------------------------------------------------------
// unitTest [-nobuf] [-mpsDir=V1] [-netlibDir=V2] [-testOsiSolverInterface]
// 
// where:
//   -nobuf: remove buffering on cout (stdout); useful to keep cout and cerr
//	 messages synchronised when redirecting output to a file or pipe.
//   -mpsDir: directory containing mps test files
//       Default value V1="../../Data/Sample"    
//   -netlibDir: directory containing netlib files
//       Default value V2="../../Data/Netlib"
//   -testOsiSolverInterface
//       If specified, then OsiSolveInterface::unitTest
//       is skipped over and not run.
//
// All parameters are optional.
//----------------------------------------------------------------

int main (int argc, const char *argv[])

{ int totalErrCnt = 0;

/*
  Start off with various bits of initialisation that don't really belong
  anywhere else.

  First off, synchronise C++ stream i/o with C stdio. This makes debugging
  output a bit more comprehensible. It still suffers from interleave of cout
  (stdout) and cerr (stderr), but -nobuf deals with that.
*/
  std::ios::sync_with_stdio() ;
/*
  Suppress an popup window that Windows shows in response to a crash. See
  note at head of file.
*/
  WindowsErrorPopupBlocker();

/*
  Process command line parameters.
*/
  std::map<std::string,std::string> parms ;

  if (processParameters(argc,argv,parms) == false)
  { return (1) ; }

  std::string mpsDir = parms["-mpsDir"] ;
  std::string netlibDir = parms["-netlibDir"] ;

try {
/*
  Test Osi{Row,Col}Cut routines.
*/
  {
    OsiCbcSolverInterface cbcSi;
    testingMessage( "Testing OsiRowCut with OsiCbcSolverInterface\n" );
    OsiRowCutUnitTest(&cbcSi,mpsDir);
  }
  {
    OsiCbcSolverInterface cbcSi;
    testingMessage( "Testing OsiColCut with OsiCbcSolverInterface\n" );
    OsiColCutUnitTest(&cbcSi,mpsDir);
  }
  {
    OsiCbcSolverInterface cbcSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiCbcSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&cbcSi,mpsDir);
  }

/*
  Run the OsiXXX class test. It's up to the OsiCbc implementor
  to decide whether or not to run OsiSolverInterfaceCommonUnitTest. Arguably
  this should be required.
*/
  testingMessage( "Testing OsiCbcSolverInterface\n" );
  OsiCbcSolverInterfaceUnitTest(mpsDir,netlibDir);

/*
  We have run the specialised unit test. Check now to see if we need to
  run through the Netlib problems.
*/
  if (parms.find("-testOsiSolverInterface") != parms.end())
  {
    // Create vector of solver interfaces
    std::vector<OsiSolverInterface*> vecSi(1, new OsiCbcSolverInterface);

    testingMessage( "Testing OsiSolverInterface on Netlib problems.\n" );
    OsiSolverInterfaceMpsUnitTest(vecSi,netlibDir);

    delete vecSi[0];
  }
  else {
    testingMessage( "***Skipped Testing of OsiCbcSolverInterface on Netlib problems***\n" );
    testingMessage( "***use -testOsiSolverInterface to run them.***\n" );
  }
} catch (CoinError& error) {
  std::cout.flush();
  std::cerr << "Caught CoinError exception: ";
  error.print(true);
  totalErrCnt += 1;
}
/*
  We're done. Report on the results.
*/
  if (totalErrCnt)
  { std::cout.flush() ;
    std::cerr
      << "Tests completed with " << totalErrCnt << " errors." << std::endl ; 
  } else
  { testingMessage("All tests completed successfully\n") ; }
  return totalErrCnt;
}
