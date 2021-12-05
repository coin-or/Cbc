// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "CoinTime.hpp"
#include "CoinFileIO.hpp"
#include "CbcModel.hpp"
#include "CbcHeuristic.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcBranchCut.hpp"
#include "CglProbing.hpp"
#include "OsiClpSolverInterface.hpp"
#include "ClpFactorization.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CbcSolver.hpp"
//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

//#############################################################################

// Display message on stdout and stderr
void testingMessage(const char *const msg)
{
  std::cout << msg;
  //cout <<endl <<"*****************************************"
  //     <<endl <<msg <<endl;
}

//#############################################################################

static inline bool CbcTestFile(const std::string name)
{
  FILE *fp = fopen(name.c_str(), "r");
  if (fp) {
    fclose(fp);
    return true;
  }
  return false;
}

//#############################################################################

bool CbcTestMpsFile(std::string &fname)
{
  if (CbcTestFile(fname)) {
    return true;
  }
  if (CbcTestFile(fname + ".mps")) {
    fname += ".mps";
    return true;
  }
  if (CbcTestFile(fname + ".MPS")) {
    fname += ".MPS";
    return true;
  }
  if (CoinFileInput::haveGzipSupport()) {
    if (CbcTestFile(fname + ".gz")) {
      return true;
    }
    if (CbcTestFile(fname + ".mps.gz")) {
      fname += ".mps";
      return true;
    }
    if (CbcTestFile(fname + ".MPS.gz")) {
      fname += ".MPS";
      return true;
    }
    if (CbcTestFile(fname + ".MPS.GZ")) {
      fname += ".MPS";
      return true;
    }
  }
  if (CoinFileInput::haveBzip2Support()) {
    if (CbcTestFile(fname + ".bz2")) {
      return true;
    }
    if (CbcTestFile(fname + ".mps.bz2")) {
      fname += ".mps";
      return true;
    }
    if (CbcTestFile(fname + ".MPS.bz2")) {
      fname += ".MPS";
      return true;
    }
    if (CbcTestFile(fname + ".MPS.BZ2")) {
      fname += ".MPS";
      return true;
    }
  }
  return false;
}
//#############################################################################
/*
  jjf: testSwitch -2 unitTest, -1 normal (==2)

  MiplibTest might be more appropriate.

  TestSwitch and stuff[6] together control how much of miplib is executed:
    For testSwitch set to:
      -3: solve 0,1 using CbcMain1
      -2: solve p0033 and p0201 only (the unit test)
      -1: solve miplib sets #0 and #1
       0: solve nothing
       k: execute sets j:k, where j is determined by the value of stuff[6]
  The last parameter of PUSH_MPS specifies the test set membership.

  For -miplib, -extra2 sets testSwitch, -extra3 sets stuff[6]. The command
    cbc -extra2 -2 -miplib
  will execute the unit test on the miplib directory.

  dirMiplib should end in the directory separator character for the platform.

  If you want to activate the row cut debugger for a given problem, change the
  last parameter of the PUSH_MPS macro for the problem to true.

  Returns 0 if all goes well, -1 if the Miplib directory is missing, otherwise
     100*(number with bad objective)+(number that exceeded node limit)
*/
int CbcClpUnitTest(const CbcModel &saveModel, const std::string &dirMiplibIn,
		   int testSwitch, const double *stuff, int argc,
		   const char ** argv,
		   int callBack(CbcModel *currentSolver, int whereFrom),
		   CbcSolverUsefulData &parameterData)
{
  // Stop Windows popup
  WindowsErrorPopupBlocker();
  unsigned int m;

  // which levels to do
  char doThisSet[40]={0};
  // default
  if (testSwitch < 0) {
    doThisSet[0]=1;
    doThisSet[1]=1;
  }
  // See if we want to use CbcMain0/1
  bool oldStyle=true;
  std::string dirMiplib = dirMiplibIn;
  // see if level indicated
  size_t found = dirMiplib.find("++level");
  size_t length = dirMiplib.size();
  int hiSet = 0;
  if (found<length) {
    // yes
    doThisSet[0]=0; doThisSet[1]=0;
    oldStyle = false;
    std::string sub = dirMiplib.substr(found+7,length-(found+8));
    dirMiplib = dirMiplib.substr(0,found+1);
    dirMiplib[found] = dirMiplibIn[length-1];
    found  = 0;
    length = sub.size();
    int level;
    while (found<length) {
      size_t foundm = sub.find("-",found);
      size_t foundc = sub.find(",",found);
      if (foundm < foundc) {
	int lo;
	std::istringstream(sub.substr(found,foundm)) >> lo;
	if (lo>30 || lo <0) {
	  std::cout << "Unable to convert to integer " << sub.substr(found)
		    << std::endl;
	  break;
	}
	found = foundm+1;
	foundc = sub.find(",",found);
	int hi;
	std::istringstream(sub.substr(found,foundc)) >> hi;
	if (hi>30 || hi <0) {
	  std::cout << "Unable to convert to integer " << sub.substr(found)
		    << std::endl;
	  break;
	}
	if (lo > hi) {
	  int i = lo;
	  lo = hi;
	  hi = i;
	}
	for (int i=lo;i<=hi;i++)
	  doThisSet[i] = 1;
      } else if (foundm > foundc) {
	std::istringstream(sub.substr(found,foundc)) >> level;
	if (level>30 || level <0) {
	  std::cout << "Unable to convert to integer " << sub.substr(found)
		    << std::endl;
	  break;
	} else {
	  doThisSet[level]=1;
	}
      } else {
	// end
	std::istringstream(sub.substr(found)) >> level;
	if (level>30 || level <0) {
	  std::cout << "Unable to convert to integer " << sub.substr(found)
		    << std::endl;
	  break;
	} else {
	  doThisSet[level]=1;
	}
      }
      if (foundc > length)
	break;
      found = foundc+1;
    }
    testSwitch = 100000;
    std::cout << " testing sets ";
    for (int i=0;i<30;i++) {
      if (doThisSet[i]) {
	std::cout << i << " ";
	if (i > 20)
	  testSwitch = 1000001; // user stuff
	hiSet = i;
      }
    }
    std::cout << std::endl;
  }
  // Do an existence check.
  std::string test1 = dirMiplib + "p0033";
  bool doTest = CbcTestMpsFile(test1);
  if (!doTest) {
    if (testSwitch >=1000000) {
      // miplib2010 or user
      test1 = dirMiplib + ((testSwitch==1000000) ? "mzzv11" : "usertest1");
      doTest = CbcTestMpsFile(test1);
      if (!doTest) {
	std::cout
	  << "Not doing miplib run as can't find mps files." << std::endl
	  << "Perhaps you're trying to read gzipped (.gz) files without libz?"
	  << std::endl;
	return (0);
      }
    }
  }
  // See if we want to use CbcMain0/1
  if (testSwitch==-3) {
    oldStyle = false;
    testSwitch=-1;
  }
  if (testSwitch >=1000000) {
    oldStyle = false;
  }
  int dfltPrecision = static_cast< int >(std::cout.precision());
  /*
  Set the range of problems to be tested. testSwitch = -2 is special and is
  picked up below.
  */
  /*
  Vectors to hold test problem names and characteristics.
*/
  std::vector< std::string > mpsName;
  std::vector< int > nRows;
  std::vector< int > nCols;
  std::vector< double > objValueC;
  std::vector< double > objValue;
  std::vector< int > testSet;
  std::vector< bool > rowCutDebugger;
/*
  A macro to make the vector creation marginally readable. Parameters are
  name, rows, columns, integer objective, continuous objective, set ID,
  row cut debugger

  To enable the row cut debugger for a given problem, change the last
  parameter to true. Don't forget to turn it off before committing changes!
*/
#define PUSH_MPS(zz_mpsName_zz,                              \
  zz_nRows_zz, zz_nCols_zz, zz_objValue_zz, zz_objValueC_zz, \
  zz_testSet_zz, zz_rcDbg_zz)                                \
  mpsName.push_back(zz_mpsName_zz);                          \
  nRows.push_back(zz_nRows_zz);                              \
  nCols.push_back(zz_nCols_zz);                              \
  objValueC.push_back(zz_objValueC_zz);                      \
  testSet.push_back(zz_testSet_zz);                          \
  objValue.push_back(zz_objValue_zz);                        \
  rowCutDebugger.push_back(zz_rcDbg_zz);
  /*
  Push the miplib problems. Except for -2 (unitTest), push all, even if we're
  not going to do all of them.
*/
  if (testSwitch == -2) {
    PUSH_MPS("p0033", 16, 33, 3089, 2520.57, 0, false);
    PUSH_MPS("p0201", 133, 201, 7615, 6875.0, 0, false);
    // PUSH_MPS("flugpl", 18, 18, 1201500, 1167185.7, 0, false);
  } else {
/*
  Load up the problem vector. Note that the row counts here include the
  objective function.
*/
    PUSH_MPS("10teams", 230, 2025, 924, 917, 1, false);
    PUSH_MPS("air03", 124, 10757, 340160, 338864.25, 0, false);
    PUSH_MPS("air04", 823, 8904, 56137, 55535.436, 2, false);
    PUSH_MPS("air05", 426, 7195, 26374, 25877.609, 2, false);
    PUSH_MPS("arki001", 1048, 1388, 7580813.0459, 7579599.80787, 7, false);
    PUSH_MPS("bell3a", 123, 133, 878430.32, 862578.64, 0, false);
    PUSH_MPS("bell5", 91, 104, 8966406.49, 8608417.95, 1, false);
    PUSH_MPS("blend2", 274, 353, 7.598985, 6.9156751140, 0, false);
    PUSH_MPS("cap6000", 2176, 6000, -2451377, -2451537.325, 1, false);
    PUSH_MPS("dano3mip", 3202, 13873, 728.1111, 576.23162474, 7, false);
    PUSH_MPS("danoint", 664, 521, 65.666667, 62.637280418, 6, false);
    PUSH_MPS("dcmulti", 290, 548, 188182, 183975.5397, 0, false);
    PUSH_MPS("dsbmip", 1182, 1886, -305.19817501, -305.19817501, 0, false);
    PUSH_MPS("egout", 98, 141, 568.101, 149.589, 0, false);
    PUSH_MPS("enigma", 21, 100, 0.0, 0.0, 0, false);
    PUSH_MPS("fast0507", 507, 63009, 174, 172.14556668, 5, false);
    PUSH_MPS("fiber", 363, 1298, 405935.18000, 156082.51759, 0, false);
    PUSH_MPS("fixnet6", 478, 878, 3983, 1200.88, 1, false);
    PUSH_MPS("flugpl", 18, 18, 1201500, 1167185.7, 0, false);
    PUSH_MPS("gen", 780, 870, 112313, 112130.0, 0, false);
    PUSH_MPS("gesa2", 1392, 1224, 25779856.372, 25476489.678, 1, false);
    PUSH_MPS("gesa2_o", 1248, 1224, 25779856.372, 25476489.678, 1, false);
    PUSH_MPS("gesa3", 1368, 1152, 27991042.648, 27833632.451, 0, false);
    PUSH_MPS("gesa3_o", 1224, 1152, 27991042.648, 27833632.451, 0, false);
    PUSH_MPS("gt2", 29, 188, 21166.000, 13460.233074, 0, false);
    PUSH_MPS("harp2", 112, 2993, -73899798.00, -74353341.502, 6, false);
    PUSH_MPS("khb05250", 101, 1350, 106940226, 95919464.0, 0, false);
    PUSH_MPS("l152lav", 97, 1989, 4722, 4656.36, 1, false);
    PUSH_MPS("lseu", 28, 89, 1120, 834.68, 0, false);
    PUSH_MPS("mas74", 13, 151, 11801.18573, 10482.79528, 3, false);
    PUSH_MPS("mas76", 12, 151, 40005.05414, 38893.9036, 2, false);
    PUSH_MPS("misc03", 96, 160, 3360, 1910., 0, false);
    PUSH_MPS("misc06", 820, 1808, 12850.8607, 12841.6, 0, false);
    PUSH_MPS("misc07", 212, 260, 2810, 1415.0, 1, false);
    PUSH_MPS("mitre", 2054, 10724, 115155, 114740.5184, 1, false);
    PUSH_MPS("mkc", 3411, 5325, -563.84601, -611.85, 7, false); 
    PUSH_MPS("mod008", 6, 319, 307, 290.9, 0, false);
    PUSH_MPS("mod010", 146, 2655, 6548, 6532.08, 0, false);
    PUSH_MPS("mod011", 4480, 10958, -54558535, -62121982.55, 2, false);
    PUSH_MPS("modglob", 291, 422, 20740508, 20430947., 2, false);
    PUSH_MPS("noswot", 182, 128, -41, -43.0, 6, false);
    PUSH_MPS("nw04", 36, 87482, 16862, 16310.66667, 1, false);
    PUSH_MPS("p0033", 16, 33, 3089, 2520.57, 0, false);
    PUSH_MPS("p0201", 133, 201, 7615, 6875.0, 0, false);
    PUSH_MPS("p0282", 241, 282, 258411, 176867.50, 0, false);
    PUSH_MPS("p0548", 176, 548, 8691, 315.29, 0, false);
    PUSH_MPS("p2756", 755, 2756, 3124, 2688.75, 0, false);
    PUSH_MPS("pk1", 45, 86, 11.0, 0.0, 2, false);
    PUSH_MPS("pp08a", 136, 240, 7350.0, 2748.3452381, 1, false);
    PUSH_MPS("pp08aCUTS", 246, 240, 7350.0, 5480.6061563, 1, false);
    PUSH_MPS("qiu", 1192, 840, -132.873137, -931.638857, 3, false);
    PUSH_MPS("qnet1", 503, 1541, 16029.692681, 14274.102667, 0, false);
    PUSH_MPS("qnet1_o", 456, 1541, 16029.692681, 12095.571667, 0, false);
    PUSH_MPS("rentacar", 6803, 9557, 30356761, 28806137.644, 0, false);
    PUSH_MPS("rgn", 24, 180, 82.1999, 48.7999, 0, false);
    PUSH_MPS("rout", 291, 556, 1077.56, 981.86428571, 3, false);
    PUSH_MPS("set1ch", 492, 712, 54537.75, 32007.73, 5, false);
    PUSH_MPS("seymour", 4944, 1372, 423, 403.84647413, 7, false);
    PUSH_MPS("stein27", 118, 27, 18, 13.0, 0, false);
    PUSH_MPS("stein45", 331, 45, 30, 22.0, 1, false);
    PUSH_MPS("swath", 884, 6805, 497.603, 334.4968581, 7, false);
    PUSH_MPS("vpm1", 234, 378, 20, 15.4167, 0, false);
    PUSH_MPS("vpm2", 234, 378, 13.75, 9.8892645972, 0, false);
    /*
      The user can add some miplib2010 models to a unit test.
      You do this by using something like 
      -dirmiplib ../miplib2010++level10-12 .... -unittest
      This will do ones with difficulty 10, 11 and 12
      other example swould be ..++level11, ++level12,14 etc 
      Difficulty 10-13 not too bad
      14-17 harder
      18 large may run out of memory when threaded
      20 infeasible problems
      21 - 28 user problems
      If the user wants to check that usertest1 has value 12345.0 then
      -DCBC_USER_UNIT_TEST1=12345.0 in build
     */
    PUSH_MPS("30n20b8", 576, 18380, 302, 1.56641,12,false);
    PUSH_MPS("acc-tight5", 3052, 1339, 0, 0,10,false);
    PUSH_MPS("aflow40b", 1442, 2728, 1168, 1005.66,12,false);
    PUSH_MPS("air04", 823, 8904, 56137, 55535.3,10,false);
    PUSH_MPS("app1-2", 53467, 26871, -41, -264.602,12,false);
    PUSH_MPS("ash608gpia-3col", 24748, 3651,1.0e50, 2,20,false); // infeasible
    PUSH_MPS("bab5", 4964, 21600, -106411.8401, -124658,12,false);
    PUSH_MPS("beasleyC3", 1750, 2500, 754, 40.4268,12,false);
    PUSH_MPS("biella1", 1203, 7328, 3.06500578e+6, 3.06004e+06,12,false);
    PUSH_MPS("bienst2", 576, 505, 54.6, 11.7241,10,false);
    PUSH_MPS("binkar10_1", 1026, 2298, 6742.2, 6637.19,10,false);
    PUSH_MPS("bley_xl1", 175620, 5831, 190, 140,12,false);
    PUSH_MPS("bnatt350", 4923, 3150, 0, 0,12,false);
    PUSH_MPS("core2536-691", 2539, 15293, 689, 688.476,11,false);
    PUSH_MPS("cov1075", 637, 120, 20, 17.1429,11,false);
    PUSH_MPS("csched010", 351, 1758, 408, 332.423,12,false);
    PUSH_MPS("danoint", 664, 521, 65.6666667, 62.6373,12,false);
    PUSH_MPS("dfn-gwin-UUM", 158, 938, 38752, 27467.3,11,false);
    PUSH_MPS("eil33-2", 32, 4516, 934.007916, 811.279,12,false);
    PUSH_MPS("eilB101", 100, 2818, 1216.92017, 1075.25,10,false);
    PUSH_MPS("enlight13", 169, 338, 71, 0,12,false);
    PUSH_MPS("enlight14", 196, 392, 1.0e50, 0,20,false); // infeasible
    PUSH_MPS("ex9", 40962, 10404, 81, 81,10,false); // likes heavy probing
    PUSH_MPS("glass4", 396, 322, 1.2000126e+09, 8.00002e+08,12,false);
    PUSH_MPS("gmu-35-40", 424, 1205, -2.49673337e+06, -2.40694e+06,12,false);
    PUSH_MPS("iis-100-0-cov", 3831, 100, 29, 16.6667,12,false);
    PUSH_MPS("iis-bupa-cov", 4803, 345, 36, 26.4972,12,false);
    PUSH_MPS("iis-pima-cov", 7201, 768, 33, 26.6204,12,false);
    PUSH_MPS("lectsched-4-obj", 14163, 7901, 4, 0,10,false);
    PUSH_MPS("m100n500k4r1", 100, 500, -25, -25,12,false);
    PUSH_MPS("macrophage", 3164, 2260, 374, 0,11,false);
    PUSH_MPS("map18", 328818, 164547, -847, -932.783,12,false);
    PUSH_MPS("map20", 328818, 164547, -922, -998.836,12,false);
    PUSH_MPS("mcsched", 2107, 1747, 211913, 193775,11,false);
    PUSH_MPS("mik-250-1-100-1", 151, 251, -66729, -79842.4,10,false);
    PUSH_MPS("mine-166-5", 8429, 830, -5.66396e+08, -8.21764e+08,10,false);
    PUSH_MPS("mine-90-10", 6270, 900, -7.843023e+08, -8.87165e+08,12,false);
    PUSH_MPS("msc98-ip", 15850, 21143, 1.9839497e+07, 1.9521e+07,12,false);
    PUSH_MPS("mspp16", 561657, 29280, 363, 341,18,false);
    PUSH_MPS("mzzv11", 9499, 10240, -21718, -22945.2,12,false);
    PUSH_MPS("n3div36", 4484, 22120, 130800, 114333,11,false);
    PUSH_MPS("n3seq24", 6044, 119856, 52200, 52000,12,false);
    PUSH_MPS("n4-3", 1236, 3596, 8993, 4080.88,12,false);
    PUSH_MPS("neos-1109824", 28979, 1520, 378, 278,11,false);
    PUSH_MPS("neos-1337307", 5687, 2840, -202319, -203124,12,false);
    PUSH_MPS("neos-1396125", 1494, 1161, 3000.04534, 388.552,11,false);
    PUSH_MPS("neos13", 20852, 1827, -95.47481, -126.178,11,false);
    PUSH_MPS("neos-1601936", 3131, 4446, 3, 1,12,false);
    PUSH_MPS("neos18", 11402, 3312, 16, 7,11,false);
    PUSH_MPS("neos-476283", 10015, 11915, 406.363, 406.245,12,false);
    PUSH_MPS("neos-686190", 3664, 3660, 6730, 5134.81,10,false);
    PUSH_MPS("neos-849702", 1041, 1737, 0, 0,11,false);
    PUSH_MPS("neos-916792", 1909, 1474, 31.870398, 26.2036,12,false);
    PUSH_MPS("neos-934278", 11495, 23123, 260, 259.5,12,false);
    PUSH_MPS("net12", 14021, 14115, 214, 17.2495,11,false);
    PUSH_MPS("netdiversion", 119589, 129180, 242, 230.8,12,false);
    PUSH_MPS("newdano", 576, 505, 65.666667, 11.7241,12,false);
    PUSH_MPS("noswot", 182, 128, -41, -43,12,false);
    PUSH_MPS("ns1208400", 4289, 2883, 2, 0,12,false);
    PUSH_MPS("ns1688347", 4191, 2685, 27, 2,12,false);
    PUSH_MPS("ns1758913", 624166, 17956, -1454.67, -1501.18,12,false);
    PUSH_MPS("ns1766074", 182, 100, 1.0e50, 5833.8,20,false); // infeasible
    PUSH_MPS("ns1830653", 2932, 1629, 20622, 6153,12,false);
    PUSH_MPS("opm2-z7-s2", 31798, 2023, -10280, -12879.7,11,false);
    PUSH_MPS("pg5_34", 225, 2600, -14339.35, -16646.6,11,false);
    PUSH_MPS("pigeon-10", 931, 490, -9000, -10000,11,false);
    PUSH_MPS("pw-myciel4", 8164, 1059, 10, 0,11,false);
    PUSH_MPS("qiu", 1192, 840, -132.873, -931.639,10,false);
    PUSH_MPS("rail507", 509, 63019, 174, 172.146,11,false);
    PUSH_MPS("ran16x16", 288, 512, 3823, 3116.43,11,false);
    PUSH_MPS("reblock67", 2523, 670, -3.4630648e+07, -3.93399e+07,12,false);
    PUSH_MPS("rmatr100-p10", 7260, 7359, 423, 360.593,10,false);
    PUSH_MPS("rmatr100-p5", 8685, 8784, 976, 762.04,10,false);
    PUSH_MPS("rmine6", 7078, 1096, -457.186, -462.306,11,false);
    PUSH_MPS("rocII-4-11", 21738, 9234, -6.65276, -11.9372,12,false);
    PUSH_MPS("rococoC10-001000", 1293, 3117, 11460, 7515.27,11,false);
    PUSH_MPS("roll3000", 2295, 1166, 12890, 11097.1,10,false);
    PUSH_MPS("satellites1-25", 5996, 9013, -5, -20,11,false);
    PUSH_MPS("sp98ic", 825, 10894, 4.49145e+08, 4.44278e+08,12,false);
    PUSH_MPS("sp98ir", 1531, 1680, 2.19677e+08, 2.16663e+08,10,false);
    PUSH_MPS("tanglegram1", 68342, 34759, 5182, 0,12,false);
    PUSH_MPS("tanglegram2", 8980, 4714, 443, 0,10,false);
    PUSH_MPS("timtab1", 171, 397, 764772, 28694,12,false);
    PUSH_MPS("triptim1", 15706, 30055, 22.8681, 22.8681,11,false);
    PUSH_MPS("unitcal_7", 48939, 25755, 1.96356e+07, 1.93876e+07,12,false);
    PUSH_MPS("vpphard", 47280, 51471, 5, -2.94558e-09,12,false);
    PUSH_MPS("zib54-UUE", 1809, 5150, 1.0334e+07, 3.87586e+06,12,false);
#ifndef CBC_USER_UNIT_TEST1
#define CBC_USER_UNIT_TEST1 -1.0e50
#endif
    PUSH_MPS("usertest1", 0, 0, CBC_USER_UNIT_TEST1, 0, 21,false);
#ifndef CBC_USER_UNIT_TEST2
#define CBC_USER_UNIT_TEST2 -1.0e50
#endif
    PUSH_MPS("usertest2", 0, 0, CBC_USER_UNIT_TEST2, 0, 22,false);
#ifndef CBC_USER_UNIT_TEST3
#define CBC_USER_UNIT_TEST3 -1.0e50
#endif
    PUSH_MPS("usertest3", 0, 0, CBC_USER_UNIT_TEST3, 0, 23,false);
#ifndef CBC_USER_UNIT_TEST4
#define CBC_USER_UNIT_TEST4 -1.0e50
#endif
    PUSH_MPS("usertest4", 0, 0, CBC_USER_UNIT_TEST4, 0, 24,false);
#ifndef CBC_USER_UNIT_TEST5
#define CBC_USER_UNIT_TEST5 -1.0e50
#endif
    PUSH_MPS("usertest5", 0, 0, CBC_USER_UNIT_TEST5, 0, 25,false);
#ifndef CBC_USER_UNIT_TEST6
#define CBC_USER_UNIT_TEST6 -1.0e50
#endif
    PUSH_MPS("usertest6", 0, 0, CBC_USER_UNIT_TEST6, 0, 26,false);
#ifndef CBC_USER_UNIT_TEST7
#define CBC_USER_UNIT_TEST7 -1.0e50
#endif
    PUSH_MPS("usertest7", 0, 0, CBC_USER_UNIT_TEST7, 0, 27,false);
#ifndef CBC_USER_UNIT_TEST8
#define CBC_USER_UNIT_TEST8 -1.0e50
#endif
    PUSH_MPS("usertest8", 0, 0, CBC_USER_UNIT_TEST8, 0, 28,false);
  }
#undef PUSH_MPS

  /*
  Normally the problems are executed in order. Define RANDOM_ORDER below to
  randomize.

  #define RANDOM_ORDER
*/
  int which[200];
  int nLoop = static_cast< int >(mpsName.size());
  assert(nLoop <= 200);
  for (int i = 0; i < nLoop; i++)
    which[i] = i;

#ifdef RANDOM_ORDER
  unsigned int iTime = static_cast< unsigned int >(CoinGetTimeOfDay() - 1.256e9);
  std::cout << "Time (seed) " << iTime << "." << std::endl;
  double sort[100];
  CoinDrand48(true, iTime);
  for (int i = 0; i < nLoop; i++)
    sort[i] = CoinDrand48();
  CoinSort_2(sort, sort + nLoop, which);
#endif

  int problemCnt = 0;
  for (m = 0; m < mpsName.size(); m++) {
    int setID = testSet[m];
    if (doThisSet[setID])
      problemCnt++;
  }

  int numberFailures = 0;
  int numberAttempts = 0;
  int numProbSolved = 0;
  double timeTaken = 0.0;

//#define CLP_FACTORIZATION_INSTRUMENT
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double timeTakenFac = 0.0;
#endif
  /*
  Open the main loop to step through the MPS problems.
  */
  for (unsigned int mw = 0; mw < mpsName.size(); mw++) {
    m = which[mw];
    int setID = testSet[m];
    // Skip if problem is not in specified problem set(s)
    if (!doThisSet[setID])
      continue;

    numberAttempts++;
    std::cout << "  processing mps file: " << mpsName[m]
              << " (" << numberAttempts << " out of "
              << problemCnt << ")" << std::endl;
    /*
  Stage 1: Read the MPS and make sure the size of the constraint matrix
	   is correct.
*/
    CbcModel *model = NULL;
    std::string fn = dirMiplib + mpsName[m];
    if (!CbcTestMpsFile(fn)) {
      std::cout << "ERROR: Cannot find MPS file " << fn << "." << std::endl;
      continue;
    }
    // Careful! We're initialising for the benefit of other code.
    CoinDrand48(true, 123456);
    double startTime = CoinCpuTime() + CoinCpuTimeJustChildren();
    if (oldStyle) {
      model = new CbcModel(saveModel);
      model->solver()->readMps(fn.c_str(), "");
    } else {
      OsiClpSolverInterface solver1;
      const char * newArgv[200];
      char replace[100];
      int newArgc = 2;
      newArgv[0] = "unitTestCbc";
      newArgv[1] = fn.c_str();
      for (int i = 3;i < argc-1; i++) {
	if (!strstr(argv[i],"++")) {
	  if (testSwitch >=1000000) {
	    // take out dextra3
	    if (strstr(argv[i],"dextra3")) {
	      i++;
	      continue;
	    }
	  }
	  newArgv[newArgc++] = argv[i];
	} else {
	  int n = strstr(argv[i],"++")-argv[i];
	  strncpy(replace,argv[i],n);
	  const char * mipname = mpsName[m].c_str();
	  int n1 = n;
	  for (int j=0;j<strlen(mipname);j++)
	    replace[n++]=mipname[j];
	  for (int j=n1+2;j<strlen(argv[i]);j++)
	    replace[n++]=argv[i][j];
	  replace[n] = '\0';
	  newArgv[newArgc++] = replace;
	  printf("Replacing %s by %s\n",argv[i],replace);
	}
      }
      /*
	Activate the row cut debugger, if requested.
      */
      if (rowCutDebugger[m] == true) {
	newArgv[newArgc++]= "-debug";
	newArgv[newArgc++]= "unitTest";
      }
      newArgv[newArgc++] = "solve";
      model = new CbcModel(solver1);
      CbcMain0(*model,parameterData);
      CbcMain1(newArgc, newArgv, *model, callBack, parameterData);
    }
    if ((model->getNumRows() != nRows[m] ||
	 model->getNumCols() != nCols[m]) && model->getNumRows())
      printf("WARNING - model has %d row, %d columns - expected %d, %d\n",
	     model->getNumRows(),model->getNumCols(),
	     nRows[m],nCols[m]);

    if (oldStyle) {

      // Higher limits for the serious problems.
      int testMaximumNodes = 200000;
      if (hiSet > 1)
	testMaximumNodes = 20000000;
      if (model->getMaximumNodes() > testMaximumNodes) {
	model->setMaximumNodes(testMaximumNodes);
      }
      /*
	Stage 2: Call solver to solve the problem.
      */
      
#ifdef CLP_FACTORIZATION_INSTRUMENT
      extern double factorization_instrument(int type);
      double facTime1 = factorization_instrument(0);
      std::cout
	<< "Factorization - initial solve " << facTime1 << " seconds."
	<< std::endl;
      timeTakenFac += facTime1;
#endif
      
      
      // Setup specific to clp
      OsiClpSolverInterface *siClp = dynamic_cast< OsiClpSolverInterface * >(model->solver());
      ClpSimplex *modelC = NULL;
      if (siClp) {
	modelC = siClp->getModelPtr();
	ClpMatrixBase *matrix = modelC->clpMatrix();
	ClpPackedMatrix *clpMatrix = dynamic_cast< ClpPackedMatrix * >(matrix);
	if (stuff && stuff[9] && clpMatrix) {
	  // vector matrix!
	  clpMatrix->makeSpecialColumnCopy();
	}
	
#ifdef JJF_ZERO
	if (clpMatrix) {
	  int numberRows = clpMatrix->getNumRows();
	  int numberColumns = clpMatrix->getNumCols();
	  double *elements = clpMatrix->getMutableElements();
	  const int *row = clpMatrix->getIndices();
	  const CoinBigIndex *columnStart = clpMatrix->getVectorStarts();
	  const int *columnLength = clpMatrix->getVectorLengths();
	  double *smallest = new double[numberRows];
	  double *largest = new double[numberRows];
	  char *flag = new char[numberRows];
	  CoinZeroN(flag, numberRows);
	  for (int i = 0; i < numberRows; i++) {
	    smallest[i] = COIN_DBL_MAX;
	    largest[i] = 0.0;
	  }
	  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	    bool isInteger = modelC->isInteger(iColumn);
	    CoinBigIndex j;
	    for (j = columnStart[iColumn];
		 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
	      int iRow = row[j];
	      double value = fabs(elements[j]);
	      if (!isInteger)
		flag[iRow] = 1;
	      smallest[iRow] = CoinMin(smallest[iRow], value);
	      largest[iRow] = CoinMax(largest[iRow], value);
	    }
	  }
	  double *rowLower = modelC->rowLower();
	  double *rowUpper = modelC->rowUpper();
	  bool changed = false;
	  for (int i = 0; i < numberRows; i++) {
	    if (flag[i] && smallest[i] > 10.0 && false) {
	      smallest[i] = 1.0 / smallest[i];
	      if (rowLower[i] > -1.0e20)
		rowLower[i] *= smallest[i];
	      if (rowUpper[i] < 1.0e20)
		rowUpper[i] *= smallest[i];
	      changed = true;
	    } else {
	      smallest[i] = 0.0;
	    }
	  }
	  if (changed) {
	    printf("SCALED\n");
	    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
	      CoinBigIndex j;
	      for (j = columnStart[iColumn];
		   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
		int iRow = row[j];
		if (smallest[iRow])
		  elements[j] *= smallest[iRow];
	      }
	    }
	  }
	  delete[] smallest;
	  delete[] largest;
	  delete[] flag;
	}
#endif // JJF_ZERO
	
	model->checkModel();
	OsiClpSolverInterface *clpSolver = dynamic_cast< OsiClpSolverInterface * >(model);
	if (clpSolver) {
	  ClpSimplex *clps = clpSolver->getModelPtr();
	  if (clps)
	    clps->setPerturbation(50);
	}
	
	modelC->tightenPrimalBounds(0.0, 0, true);
	model->initialSolve();
	if (modelC->dualBound() == 1.0e10) {
	  // user did not set - so modify
	  // get largest scaled away from bound
	  ClpSimplex temp = *modelC;
	  temp.dual(0, 7);
	  double largestScaled = 1.0e-12;
	  double largest = 1.0e-12;
	  int numberRows = temp.numberRows();
	  const double *rowPrimal = temp.primalRowSolution();
	  const double *rowLower = temp.rowLower();
	  const double *rowUpper = temp.rowUpper();
	  const double *rowScale = temp.rowScale();
	  int iRow;
	  for (iRow = 0; iRow < numberRows; iRow++) {
	    double value = rowPrimal[iRow];
	    double above = value - rowLower[iRow];
	    double below = rowUpper[iRow] - value;
	    if (above < 1.0e12) {
	      largest = CoinMax(largest, above);
	    }
	    if (below < 1.0e12) {
	      largest = CoinMax(largest, below);
	    }
	    if (rowScale) {
	      double multiplier = rowScale[iRow];
	      above *= multiplier;
	      below *= multiplier;
	    }
	    if (above < 1.0e12) {
	      largestScaled = CoinMax(largestScaled, above);
	    }
	    if (below < 1.0e12) {
	      largestScaled = CoinMax(largestScaled, below);
	    }
	  }
	  
	  int numberColumns = temp.numberColumns();
	  const double *columnPrimal = temp.primalColumnSolution();
	  const double *columnLower = temp.columnLower();
	  const double *columnUpper = temp.columnUpper();
	  const double *columnScale = temp.columnScale();
	  int iColumn;
	  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
	    double value = columnPrimal[iColumn];
	    double above = value - columnLower[iColumn];
	    double below = columnUpper[iColumn] - value;
	    if (above < 1.0e12) {
	      largest = CoinMax(largest, above);
	    }
	    if (below < 1.0e12) {
	      largest = CoinMax(largest, below);
	    }
	    if (columnScale) {
	      double multiplier = 1.0 / columnScale[iColumn];
	      above *= multiplier;
	      below *= multiplier;
	    }
	    if (above < 1.0e12) {
	      largestScaled = CoinMax(largestScaled, above);
	    }
	    if (below < 1.0e12) {
	      largestScaled = CoinMax(largestScaled, below);
	    }
	  }
	  std::cout << "Largest (scaled) away from bound " << largestScaled
		    << " unscaled " << largest << std::endl;
#ifdef JJF_ZERO
	  modelC->setDualBound(CoinMax(1.0001e8,
				       CoinMin(1000.0 * largestScaled, 1.00001e10)));
#else
	  modelC->setDualBound(CoinMax(1.0001e9,
				       CoinMin(1000.0 * largestScaled, 1.0001e10)));
#endif
	}
      } // end clp-specific setup
      /*
	Cut passes: For small models (n < 500) always do 100 passes, if possible
	(-100). For larger models, use minimum drop to stop (100, 20).
      */
      model->setMinimumDrop(CoinMin(5.0e-2,
				    fabs(model->getMinimizationObjValue()) * 1.0e-3 + 1.0e-4));
      if (CoinAbs(model->getMaximumCutPassesAtRoot()) <= 100) {
	if (model->getNumCols() < 500) {
	  model->setMaximumCutPassesAtRoot(-100);
	} else if (model->getNumCols() < 5000) {
	  model->setMaximumCutPassesAtRoot(100);
	} else {
	  model->setMaximumCutPassesAtRoot(20);
	}
      }
      // If defaults then increase trust for small models
      if (model->numberStrong() == 5 && model->numberBeforeTrust() == 10) {
	int numberColumns = model->getNumCols();
	if (numberColumns <= 50) {
	  model->setNumberBeforeTrust(1000);
	} else if (numberColumns <= 100) {
	  model->setNumberBeforeTrust(100);
	} else if (numberColumns <= 300) {
	  model->setNumberBeforeTrust(50);
	}
      }
      //if (model->getNumCols()>=500) {
      // switch off Clp stuff
      //model->setFastNodeDepth(-1);
      //}
    /*
      Activate the row cut debugger, if requested.
    */
      if (rowCutDebugger[m] == true) {
	std::string probName;
	model->solver()->getStrParam(OsiProbName, probName);
	model->solver()->activateRowCutDebugger(probName.c_str());
	if (model->solver()->getRowCutDebugger())
	  std::cout << "Row cut debugger activated for ";
	else
	  std::cout << "Failed to activate row cut debugger for ";
	std::cout << mpsName[m] << "." << std::endl;
      }
      setCutAndHeuristicOptions(*model);
      /*
	More clp-specific setup.
      */
      if (siClp) {
#ifdef CLP_MULTIPLE_FACTORIZATIONS
	if (!modelC->factorization()->isDenseOrSmall()) {
	  int denseCode = stuff ? static_cast< int >(stuff[4]) : -1;
	  int smallCode = stuff ? static_cast< int >(stuff[10]) : -1;
	  if (stuff && stuff[8] >= 1) {
	    if (denseCode < 0)
	      denseCode = 40;
	    if (smallCode < 0)
	      smallCode = 40;
	  }
	  if (denseCode > 0)
	    modelC->factorization()->setGoDenseThreshold(denseCode);
	  if (smallCode > 0)
	    modelC->factorization()->setGoSmallThreshold(smallCode);
	  if (denseCode >= modelC->numberRows()) {
	    //printf("problem going dense\n");
	    //modelC->factorization()->goDenseOrSmall(modelC->numberRows());
	  }
	}
#endif
	if (stuff && stuff[8] >= 1) {
	  printf("Fast node size Columns %d rows %d - depth %d\n",
		 modelC->numberColumns(), modelC->numberRows(),
		 model->fastNodeDepth());
	  if (modelC->numberColumns() + modelC->numberRows() <= 10000 && model->fastNodeDepth() == -1)
	    model->setFastNodeDepth(-10 /*-9*/);
	}
      }
#ifdef CONFLICT_CUTS
      {
	model->setCutoffAsConstraint(true); // very slow on bell5 ??
	int moreOptions = model->moreSpecialOptions();
	model->setMoreSpecialOptions(moreOptions | 4194304);
      }
#endif
      /*
	Finally, the actual call to solve the MIP with branch-and-cut.
      */
      model->setDblParam( CbcModel::CbcAllowableGap, 0.0001 );
      model->branchAndBound();
      
#ifdef CLP_FACTORIZATION_INSTRUMENT
      double facTime = factorization_instrument(0);
      std::cout << "Factorization " << facTime << " seconds." << std::endl,
	timeTakenFac += facTime;
#endif
    }
    /*
      Stage 3: Do the statistics and check the answer.
    */
    double timeOfSolution = CoinCpuTime() + CoinCpuTimeJustChildren() - startTime;
    std::cout
      << "Cuts at root node changed objective from "
      << model->getContinuousObjective() << " to "
      << model->rootObjectiveAfterCuts() << std::endl;
    int numberGenerators = model->numberCutGenerators();
    for (int iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
      CbcCutGenerator *generator = model->cutGenerator(iGenerator);
#ifdef CLIQUE_ANALYSIS
#ifndef CLP_INVESTIGATE
      CglImplication *implication = dynamic_cast< CglImplication * >(generator->generator());
      if (implication)
        continue;
#endif
#endif
      std::cout
        << generator->cutGeneratorName() << " was tried "
        << generator->numberTimesEntered() << " times and created "
        << generator->numberCutsInTotal() << " cuts of which "
        << generator->numberCutsActive()
        << " were active after adding rounds of cuts";
      if (generator->timing())
        std::cout << " (" << generator->timeInCutGenerator() << " seconds)";
      std::cout << "." << std::endl;
    }
    std::cout
      << model->getNumberHeuristicSolutions()
      << " solutions found by heuristics." << std::endl;
    int numberHeuristics = model->numberHeuristics();
    for (int iHeuristic = 0; iHeuristic < numberHeuristics; iHeuristic++) {
      CbcHeuristic *heuristic = model->heuristic(iHeuristic);
      if (heuristic->numRuns()) {
        std::cout
          << heuristic->heuristicName() << " was tried "
          << heuristic->numRuns() << " times out of "
          << heuristic->numCouldRun() << " and created "
          << heuristic->numberSolutionsFound() << " solutions." << std::endl;
      }
    }
    /*
  Check for the correct answer.
*/
    double objActual = model->getObjValue();
    double objExpect = objValue[m];
    double tolerance = CoinMin(fabs(objActual), fabs(objExpect));
    tolerance = CoinMax(1.0e-4, 1.0e-5 * tolerance);
    if (!model->status()) {

      //CoinRelFltEq eq(1.0e-3) ;

      std::cout
        << "cbc_clp (" << mpsName[m] << ") "
        << std::setprecision(10) << objActual;
      if (fabs(objActual - objExpect) < tolerance) {
        std::cout << std::setprecision(dfltPrecision) << "; okay";
        numProbSolved++;
      } else if (objExpect!=-1.0e50) {
        std::cout
          << " != " << objExpect << std::setprecision(dfltPrecision)
          << "; error = " << fabs(objExpect - objActual);
        numberFailures++;
        //#ifdef COIN_DEVELOP
        //abort();
        //#endif
      } else {
        std::cout
          << " - user model ";
      }
    } else {
      std::cout
        << "cbc_clp (" << mpsName[m] << ") status not optimal; "
        << "assuming too many nodes";
      if (fabs(objActual - objExpect) < tolerance) {
	std::cout
	  << " (on the bright side solution is correct) ";
      }
    }
    timeTaken += timeOfSolution;
    std::cout
      << " -- (" << model->getNodeCount() << " n / "
      << model->getIterationCount() << " i / "
      << timeOfSolution << " s) (subtotal " << timeTaken << " seconds)"
      << std::endl << std::flush;
    delete model;
  }
  /*
  End main loop on MPS problems. Print a summary and calculate the return
  value.
*/
  int returnCode = 0;
  std::cout
    << "cbc_clp solved " << numProbSolved << " out of " << numberAttempts;
  int numberOnNodes = numberAttempts - numProbSolved - numberFailures;
  if (numberFailures || numberOnNodes) {
    if (numberOnNodes) {
      std::cout << " (" << numberOnNodes << " stopped on nodes)";
      returnCode = numberOnNodes;
    }
    if (numberFailures) {
      std::cout << " (" << numberFailures << " gave bad answer!)";
      returnCode += 100 * numberFailures;
    }
  }
  std::cout
    << " and took " << timeTaken << " seconds." << std::endl;

  if (testSwitch == -2) {
    if (numberFailures || numberOnNodes) {
      std::cout << "****** Unit Test failed." << std::endl;
      std::cerr << "****** Unit Test failed." << std::endl;
    } else {
      std::cerr << "****** Unit Test succeeded." << std::endl;
    }
  }
#ifdef CLP_FACTORIZATION_INSTRUMENT
  std::cout
    << "Total factorization time " << timeTakenFac << "seconds." << std::endl;
#endif
  return (returnCode);
}
