// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <deque>

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
  -dirMiplib <dir>++levelA-B,C chooses the sets. In a MIPLIB 2017 benchmark
  directory the sets are the levels of the MIPLIB 2017 table below, and the
  default is levels 1 and 2.

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
		   int testSwitch, const double *stuff,
		   std::deque<std::string> originalInputQueue,
		   int callBack(CbcModel *currentSolver, int whereFrom),
		   CbcParameters &parameters)
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
  bool levelGiven = false;
  if (found<length) {
    // yes
    doThisSet[0]=0; doThisSet[1]=0;
    oldStyle = false;
    levelGiven = true;
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
  /*
    The MIPLIB 2017 benchmark set (miplib.zib.de, .mps.gz files): recognised
    by the directory's name or by one of its instances that isn't in MIPLIB 3
    or 2010. It has its own levels (see the table below); without ++level,
    levels 1 and 2 are done.
  */
  bool miplib2017 = false;
  if (!doTest && testSwitch < 1000000) {
    std::string lower = dirMiplib;
    for (size_t i = 0; i < lower.size(); i++)
      lower[i] = static_cast<char>(tolower(lower[i]));
    miplib2017 = lower.find("miplib2017") != std::string::npos;
    const char *only2017[] = {"neos-3083819-nubu", "neos-2987310-joes",
                              "neos-5104907-jarama"};
    for (int i = 0; i < 3 && !miplib2017; i++) {
      std::string name = dirMiplib + only2017[i];
      miplib2017 = CbcTestMpsFile(name);
    }
    if (miplib2017) {
      doTest = true;
      if (!levelGiven) {
        doThisSet[0] = 0;
        doThisSet[1] = 1;
        doThisSet[2] = 1;
      }
      std::cout << "MIPLIB 2017 benchmark set: testing levels";
      for (int i = 0; i < 30; i++) {
        if (doThisSet[i])
          std::cout << " " << i;
      }
      std::cout << std::endl;
    }
  }
  if (!doTest) {
    if (!strstr(dirMiplib.c_str(),"miplib"))
      testSwitch = 1000001; // user test set
    if (testSwitch >=1000000) {
      // miplib2010 or user
      test1 = dirMiplib + ((testSwitch==1000000) ? "mzzv11" : "testset");
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
  } else if (miplib2017) {
    /*
      The 239 MIPLIB 2017 benchmark instances with a known optimum or proven
      infeasible (supportcase22, still open, is left out). Rows and columns
      are as read, without the objective. The continuous objective is the LP
      relaxation as read (0 where it wasn't recorded); the test doesn't use it.
      Infeasible instances have 1.0e50.

      Levels, from one-thread, one-hour runs with default settings (Cbc master
      of November 2026 and Cbc 2.10, on a cluster node; times are master's):
        1  master solves in under 10 s
        2  10-30 s (1 and 2, the default, take about 3 minutes)
        3  30-120 s
        4  120-600 s
        5  600-3600 s, or only 2.10 solves within the hour
        6  neither solves within the hour
        7  infeasible
      e.g. -dirMiplib ~/miplib2017++level3-4 -unitTest
    */
    // BEGIN MIPLIB 2017 TABLE
    PUSH_MPS("30n20b8", 576, 18380, 302.0, 1.56641, 5, false);
    PUSH_MPS("50v-10", 233, 2013, 3311.1799841, 2879.07, 6, false);
    PUSH_MPS("academictimetablesmall", 23294, 28926, 0.0, 0, 6, false);
    PUSH_MPS("air05", 426, 7195, 26374.0, 25877.6, 2, false);
    PUSH_MPS("app1-1", 4926, 2480, -3.0, -22.8244, 1, false);
    PUSH_MPS("app1-2", 53467, 26871, -41.0, -264.602, 6, false);
    PUSH_MPS("assign1-5-8", 161, 156, 211.999999999998, 183.363, 6, false);
    PUSH_MPS("atlanta-ip", 21732, 48738, 90.009878614, 81.2432, 6, false);
    PUSH_MPS("b1c1s1", 3904, 3872, 24544.25, 2866.86, 6, false);
    PUSH_MPS("bab2", 17245, 147912, -357544.3115, -432166, 6, false);
    PUSH_MPS("bab6", 29904, 114240, -284248.2307000001, -321641, 6, false);
    PUSH_MPS("beasleyC3", 1750, 2500, 753.9999999999128, 40.4268, 6, false);
    PUSH_MPS("binkar10_1", 1026, 2298, 6742.1998835, 6637.19, 5, false);
    PUSH_MPS("blp-ar98", 1128, 16021, 6205.2147104, 5891.23, 6, false);
    PUSH_MPS("blp-ic98", 717, 13640, 4491.44758395, 4331.17, 6, false);
    PUSH_MPS("bnatt400", 5614, 3600, 1.0, 0, 6, false);
    PUSH_MPS("bnatt500", 7029, 4500, 1.0e50, 0, 7, false);
    PUSH_MPS("bppc4-08", 111, 1456, 53.0, 51.59, 6, false);
    PUSH_MPS("brazil3", 14646, 23968, 24.0, 1.99999, 6, false);
    PUSH_MPS("buildingenergy", 277594, 154978, 33283.853236, 33246.2, 6, false);
    PUSH_MPS("cbs-cta", 10112, 24793, 0.0, 0, 2, false);
    PUSH_MPS("chromaticindex1024-7", 67583, 73728, 4.0, 3, 6, false);
    PUSH_MPS("chromaticindex512-7", 33791, 36864, 4.0, 3, 6, false);
    PUSH_MPS("cmflsp50-24-8-8", 3520, 16392, 55789389.886, 5.46331e+07, 6, false);
    PUSH_MPS("CMS750_4", 16381, 11697, 252.0, 250, 6, false);
    PUSH_MPS("co-100", 2187, 48417, 2639942.06, 917102, 6, false);
    PUSH_MPS("cod105", 1024, 1024, -12.0, -18.2857, 6, false);
    PUSH_MPS("comp07-2idx", 21235, 17264, 6.0, 0, 5, false);
    PUSH_MPS("comp21-2idx", 14038, 10863, 74.0, 0, 6, false);
    PUSH_MPS("cost266-UUE", 1446, 4161, 25148940.55999998, 2.01615e+07, 6, false);
    PUSH_MPS("cryptanalysiskb128n5obj14", 98021, 48950, 1.0e50, 0, 7, false);
    PUSH_MPS("cryptanalysiskb128n5obj16", 98021, 48950, 0.0, 0, 6, false);
    PUSH_MPS("csched007", 351, 1758, 350.9999999999955, 269.252, 6, false);
    PUSH_MPS("csched008", 351, 1536, 173.0, 171, 6, false);
    PUSH_MPS("cvs16r128-89", 4633, 3472, -97.0, -128, 6, false);
    PUSH_MPS("dano3_3", 3202, 13873, 576.34463303, 576.232, 3, false);
    PUSH_MPS("dano3_5", 3202, 13873, 576.9249159565619, 576.232, 4, false);
    PUSH_MPS("decomp2", 10765, 14387, -160.0, -160, 1, false);
    PUSH_MPS("drayage-100-23", 4630, 11090, 103333.87407, 11714.1, 1, false);
    PUSH_MPS("drayage-25-23", 4630, 11090, 101282.647018, 11714.1, 6, false);
    PUSH_MPS("dws008-01", 6064, 11096, 37412.60458794508, 584.5, 6, false);
    PUSH_MPS("eil33-2", 32, 4516, 934.007915999999, 811.279, 3, false);
    PUSH_MPS("eilA101-2", 100, 65832, 880.9201079999999, 803.374, 6, false);
    PUSH_MPS("enlight_hard", 100, 200, 37.0, 0, 3, false);
    PUSH_MPS("ex10", 69608, 17680, 100.0, 100, 6, false);
    PUSH_MPS("ex9", 40962, 10404, 81.0, 81, 6, false);
    PUSH_MPS("exp-1-500-5-5", 550, 990, 65887.0, 28427, 6, false);
    PUSH_MPS("fast0507", 507, 63009, 174.0, 172.146, 4, false);
    PUSH_MPS("fastxgemm-n2r6s0t2", 5998, 784, 230.0, 27, 6, false);
    PUSH_MPS("fhnw-binpack4-4", 620, 520, 1.0e50, 0, 7, false);
    PUSH_MPS("fhnw-binpack4-48", 4480, 3710, 0.0, 0, 6, false);
    PUSH_MPS("fiball", 3707, 34219, 138.0, 137.069, 3, false);
    PUSH_MPS("gen-ip002", 24, 41, -4783.733392, -4840.54, 4, false);
    PUSH_MPS("gen-ip054", 27, 30, 6840.96564179, 6765.21, 4, false);
    PUSH_MPS("germanrr", 10779, 10813, 47095869.649, 4.59801e+07, 6, false);
    PUSH_MPS("gfd-schedulen180f7d50m30k18", 457985, 227535, 1.0, 1, 6, false);
    PUSH_MPS("glass-sc", 6119, 214, 23.0, 14.0803, 6, false);
    PUSH_MPS("glass4", 396, 322, 1200012599.972384, 8.00002e+08, 6, false);
    PUSH_MPS("gmu-35-40", 424, 1205, -2406733.3688, -2.40694e+06, 6, false);
    PUSH_MPS("gmu-35-50", 435, 1919, -2607958.33, -2.60807e+06, 6, false);
    PUSH_MPS("graph20-20-1rand", 5587, 2183, -9.0, -37, 6, false);
    PUSH_MPS("graphdraw-domain", 865, 254, 19685.99997550038, 12672, 6, false);
    PUSH_MPS("h80x6320d", 6558, 12640, 6382.099048245999, 5325.16, 6, false);
    PUSH_MPS("highschool1-aigio", 92568, 320404, 0.0, 2.39234e-13, 6, false);
    PUSH_MPS("hypothyroid-k1", 5195, 2602, -2851.0, -2902.85, 3, false);
    PUSH_MPS("ic97_potential", 1046, 728, 3941.99993090225, 3868, 6, false);
    PUSH_MPS("icir97_tension", 1203, 2494, 6375.0, 6302, 6, false);
    PUSH_MPS("irish-electricity", 104259, 61728, 3723497.591396, 2.45479e+06, 6, false);
    PUSH_MPS("irp", 39, 20315, 12159.49283539698, 12123.5, 2, false);
    PUSH_MPS("istanbul-no-cutoff", 20346, 5282, 204.08170701, 51.11, 5, false);
    PUSH_MPS("k1mushroom", 16419, 8211, -3288.0, -4165.06, 6, false);
    PUSH_MPS("lectsched-5-obj", 38884, 21805, 24.0, 0, 6, false);
    PUSH_MPS("leo1", 593, 6731, 404227536.16, 3.88573e+08, 6, false);
    PUSH_MPS("leo2", 593, 11100, 404077441.12, 3.86421e+08, 6, false);
    PUSH_MPS("lotsize", 1920, 2985, 1480195.0, 348385, 6, false);
    PUSH_MPS("mad", 51, 220, 0.0268, 0, 6, false);
    PUSH_MPS("map10", 328818, 164547, -495.0, -602.176, 5, false);
    PUSH_MPS("map16715-04", 328818, 164547, -111.0, -296.268, 6, false);
    PUSH_MPS("markshare2", 7, 74, 1.0, 0, 6, false);
    PUSH_MPS("markshare_4_0", 4, 34, 1.0, 0, 3, false);
    PUSH_MPS("mas74", 13, 151, 11801.18572, 10482.8, 4, false);
    PUSH_MPS("mas76", 12, 151, 40005.05398999999, 38893.9, 3, false);
    PUSH_MPS("mc11", 1920, 3040, 11688.99999999966, 608.844, 6, false);
    PUSH_MPS("mcsched", 2107, 1747, 211913.0, 193775, 5, false);
    PUSH_MPS("mik-250-20-75-4", 195, 270, -52301.0, -61651.2, 3, false);
    PUSH_MPS("milo-v12-6-r2-40-1", 5628, 2688, 326481.14282799, 205532, 6, false);
    PUSH_MPS("momentum1", 42680, 5174, 109143.4935, 72793.3, 6, false);
    PUSH_MPS("mushroom-best", 8580, 8468, 0.0553337612, 1.26829e-05, 6, false);
    PUSH_MPS("mzzv11", 9499, 10240, -21718.0, -22945.2, 4, false);
    PUSH_MPS("mzzv42z", 10460, 11717, -20540.0, -21623, 3, false);
    PUSH_MPS("n2seq36q", 2565, 22480, 52200.0, 52000, 6, false);
    PUSH_MPS("n3div36", 4484, 22120, 130800.0, 114333, 6, false);
    PUSH_MPS("n5-3", 1062, 2550, 8104.999999993999, 2883.82, 6, false);
    PUSH_MPS("neos-1122047", 57791, 5100, 161.0, 161, 1, false);
    PUSH_MPS("neos-1171448", 13206, 4914, -309.0, -309, 3, false);
    PUSH_MPS("neos-1171737", 4179, 2340, -195.0, -195, 3, false);
    PUSH_MPS("neos-1354092", 3135, 13702, 46.0, 13.2857, 6, false);
    PUSH_MPS("neos-1445765", 2147, 20617, -17783.0, -24699.2, 6, false);
    PUSH_MPS("neos-1456979", 6770, 4605, 176.0, 154, 6, false);
    PUSH_MPS("neos-1582420", 10180, 10100, 90.99999999999996, 87.5761, 2, false);
    PUSH_MPS("neos-2075418-temuka", 349602, 122304, 1.0e50, 0, 7, false);
    PUSH_MPS("neos-2657525-crna", 342, 524, 1.810748, 0, 6, false);
    PUSH_MPS("neos-2746589-doon", 31530, 50936, 2008.199999999999, 1985.75, 6, false);
    PUSH_MPS("neos-2978193-inde", 396, 20800, -2.38806169, -2.41852, 6, false);
    PUSH_MPS("neos-2987310-joes", 29015, 27837, -607702988.3, -6.10797e+08, 2, false);
    PUSH_MPS("neos-3004026-krka", 12545, 17030, 0.0, 0, 6, false);
    PUSH_MPS("neos-3024952-loue", 3705, 3255, 26756.0, 22882, 6, false);
    PUSH_MPS("neos-3046615-murg", 498, 274, 1600.0, 192, 6, false);
    PUSH_MPS("neos-3083819-nubu", 4725, 8644, 6307996.0, 6.29058e+06, 1, false);
    PUSH_MPS("neos-3216931-puriri", 5989, 3555, 71320.0, 58817.4, 6, false);
    PUSH_MPS("neos-3381206-awhea", 479, 2375, 453.0, 415.24, 5, false);
    PUSH_MPS("neos-3402294-bobin", 591076, 2904, 0.06724999999999949, 0, 4, false);
    PUSH_MPS("neos-3402454-bohle", 2897380, 2904, 1.0e50, 0, 7, false);
    PUSH_MPS("neos-3555904-turama", 146493, 37461, -34.7, -41.45, 6, false);
    PUSH_MPS("neos-3627168-kasai", 1655, 1462, 988585.6199999998, 945808, 6, false);
    PUSH_MPS("neos-3656078-kumeu", 17656, 14870, -13172.2, -18774.7, 6, false);
    PUSH_MPS("neos-3754480-nidda", 402, 253, 12941.73838561778, -1.21692e+06, 6, false);
    PUSH_MPS("neos-3988577-wolgan", 44662, 25870, 1.0e50, 119, 7, false);
    PUSH_MPS("neos-4300652-rahue", 76992, 33003, 2.1416, 0, 6, false);
    PUSH_MPS("neos-4338804-snowy", 1701, 1344, 1471.0, 1447, 6, false);
    PUSH_MPS("neos-4387871-tavua", 4554, 4004, 33.384729927, 10.1438, 6, false);
    PUSH_MPS("neos-4413714-turia", 2303, 190402, 45.3701670199998, 43.4747, 3, false);
    PUSH_MPS("neos-4532248-waihi", 167322, 86842, 61.59999999999999, 0, 6, false);
    PUSH_MPS("neos-4647030-tutaki", 8382, 12600, 27265.70599999996, 26911.6, 6, false);
    PUSH_MPS("neos-4722843-widden", 113555, 77723, 25009.662227, -1603.54, 6, false);
    PUSH_MPS("neos-4738912-atrato", 1947, 6216, 283627956.595, 6.51435e+07, 5, false);
    PUSH_MPS("neos-4763324-toguru", 106954, 53593, 1613.03884585, 1120.99, 6, false);
    PUSH_MPS("neos-4954672-berkel", 1848, 1533, 2612710.0, 1.15023e+06, 6, false);
    PUSH_MPS("neos-5049753-cuanza", 322248, 242736, 561.9999971689, 0, 6, false);
    PUSH_MPS("neos-5052403-cygnet", 38268, 32868, 182.0, 179.5, 6, false);
    PUSH_MPS("neos-5093327-huahum", 51840, 40640, 6259.999997125895, 3934, 6, false);
    PUSH_MPS("neos-5104907-jarama", 489818, 345856, 935.0, 0, 6, false);
    PUSH_MPS("neos-5107597-kakapo", 6498, 3114, 3644.99999999952, 0, 6, false);
    PUSH_MPS("neos-5114902-kasavu", 961170, 710164, 655.0, 0, 6, false);
    PUSH_MPS("neos-5188808-nattai", 29452, 14544, 0.110283622999984, 0, 6, false);
    PUSH_MPS("neos-5195221-niemur", 42256, 14546, 0.0038354326, 0, 6, false);
    PUSH_MPS("neos-631710", 169576, 167056, 203.0, 188.25, 6, false);
    PUSH_MPS("neos-662469", 1085, 18235, 184379.9999999999, 184341, 6, false);
    PUSH_MPS("neos-787933", 1897, 236376, 30.0, 3, 6, false);
    PUSH_MPS("neos-827175", 14187, 32504, 112.00152, 112.002, 1, false);
    PUSH_MPS("neos-848589", 1484, 550539, 2351.40309999697, 0, 6, false);
    PUSH_MPS("neos-860300", 850, 1385, 3200.999999999998, 1666.07, 3, false);
    PUSH_MPS("neos-873061", 93360, 175288, 113.6562385063, 5.52781, 6, false);
    PUSH_MPS("neos-911970", 107, 888, 54.76, 23.26, 6, false);
    PUSH_MPS("neos-933966", 12047, 31762, 318.0, 318, 4, false);
    PUSH_MPS("neos-950242", 34224, 5760, 4.0, 1, 5, false);
    PUSH_MPS("neos-957323", 3757, 57756, -237.7566815, -237.758, 3, false);
    PUSH_MPS("neos-960392", 4744, 59376, -238.0, -238, 3, false);
    PUSH_MPS("neos17", 486, 535, 0.1500025774, 0.000681499, 3, false);
    PUSH_MPS("neos5", 63, 63, 15.0, 13, 5, false);
    PUSH_MPS("neos8", 46324, 23228, -3719.0, -3773.51, 3, false);
    PUSH_MPS("neos859080", 164, 160, 1.0e50, 1, 7, false);
    PUSH_MPS("net12", 14021, 14115, 214.0, 17.2495, 5, false);
    PUSH_MPS("netdiversion", 119589, 129180, 242.0, 230.8, 6, false);
    PUSH_MPS("nexp-150-20-8-5", 4620, 20115, 231.0, 18.2268, 6, false);
    PUSH_MPS("ns1116954", 131991, 12648, 0.0, 0, 6, false);
    PUSH_MPS("ns1208400", 4289, 2883, 2.0, 0, 5, false);
    PUSH_MPS("ns1644855", 40698, 30200, -1524.33333333333, -1524.33, 5, false);
    PUSH_MPS("ns1760995", 615388, 17956, -549.21438505, -1379.17, 6, false);
    PUSH_MPS("ns1830653", 2932, 1629, 20622.0, 6153, 6, false);
    PUSH_MPS("ns1952667", 41, 13264, 0.0, 0, 6, false);
    PUSH_MPS("nu25-pr12", 2313, 5868, 53904.99999999999, 52757.5, 2, false);
    PUSH_MPS("nursesched-medium-hint03", 14062, 34248, 115.0, 33.4685, 6, false);
    PUSH_MPS("nursesched-sprint02", 3522, 10250, 57.99999999999999, 54.4167, 3, false);
    PUSH_MPS("nw04", 36, 87482, 16862.0, 16310.7, 1, false);
    PUSH_MPS("opm2-z10-s4", 160633, 6250, -33269.0, -46921.5, 6, false);
    PUSH_MPS("p200x1188c", 1388, 2376, 15078.0, 5678.61, 5, false);
    PUSH_MPS("peg-solitaire-a3", 4587, 4552, 1.0, 1, 6, false);
    PUSH_MPS("pg", 125, 2700, -8674.34260712, -11824.7, 3, false);
    PUSH_MPS("pg5_34", 225, 2600, -14339.35345, -16646.6, 5, false);
    PUSH_MPS("physiciansched3-3", 266227, 79555, 2623271.326667, 0, 6, false);
    PUSH_MPS("physiciansched6-2", 168336, 111827, 49324.0, 43012.5, 2, false);
    PUSH_MPS("piperout-08", 14589, 10399, 125054.9999999999, 99128.9, 4, false);
    PUSH_MPS("piperout-27", 18442, 11659, 8123.999999999973, 6970, 3, false);
    PUSH_MPS("pk1", 45, 86, 11.0, 0, 3, false);
    PUSH_MPS("proteindesign121hz512p9", 301, 159145, 1473.0, 1423.9, 6, false);
    PUSH_MPS("proteindesign122trx11p8", 254, 127326, 1747.0, 1720.46, 6, false);
    PUSH_MPS("qap10", 1820, 4150, 339.9999999983871, 332.566, 3, false);
    PUSH_MPS("radiationm18-12-05", 40935, 40623, 17566.0, 17554.5, 6, false);
    PUSH_MPS("radiationm40-10-02", 173603, 172013, 155328.0, 155307, 6, false);
    PUSH_MPS("rail01", 46843, 117527, -70.5699643, -92.0873, 6, false);
    PUSH_MPS("rail02", 95791, 270869, -200.4499077, -206.61, 6, false);
    PUSH_MPS("rail507", 509, 63019, 174.0, 172.146, 5, false);
    PUSH_MPS("ran14x18-disj-8", 447, 504, 3712.0, 3444.42, 6, false);
    PUSH_MPS("rd-rplusc-21", 125899, 622, 165395.275295, 100, 6, false);
    PUSH_MPS("reblock115", 4735, 1150, -36800603.2332, -3.9366e+07, 6, false);
    PUSH_MPS("rmatr100-p10", 7260, 7359, 423.0, 360.593, 3, false);
    PUSH_MPS("rmatr200-p5", 37617, 37816, 4521.0, 3283.65, 6, false);
    PUSH_MPS("rocI-4-11", 10883, 6839, -6020203.0, -1.11111e+07, 6, false);
    PUSH_MPS("rocII-5-11", 26897, 11523, -6.675504731538, -11.9497, 6, false);
    PUSH_MPS("rococoB10-011000", 1667, 4456, 19449.0, 8350.2, 6, false);
    PUSH_MPS("rococoC10-001000", 1293, 3117, 11460.0, 7515.27, 5, false);
    PUSH_MPS("roi2alpha3n4", 1251, 6816, -63.20849503, -83.881, 5, false);
    PUSH_MPS("roi5alpha10n8", 4665, 106150, -52.322274351, -91.9164, 6, false);
    PUSH_MPS("roll3000", 2295, 1166, 12889.999992, 11097.1, 5, false);
    PUSH_MPS("s100", 14733, 364417, -0.1697235270583, -0.171087, 6, false);
    PUSH_MPS("s250r10", 10962, 273142, -0.1717804834232, -0.172677, 5, false);
    PUSH_MPS("satellites2-40", 20916, 35378, -19.0, -30, 6, false);
    PUSH_MPS("satellites2-60-fs", 16516, 35378, -19.0000000001, -30, 6, false);
    PUSH_MPS("savsched1", 295989, 328575, 3217.7, 217.404, 6, false);
    PUSH_MPS("sct2", 2151, 5885, -230.9891623, -231.117, 6, false);
    PUSH_MPS("seymour", 4944, 1372, 423.0, 403.846, 6, false);
    PUSH_MPS("seymour1", 4944, 1372, 410.76370139, 403.846, 4, false);
    PUSH_MPS("sing326", 50781, 55156, 7753674.85376, 7.71767e+06, 6, false);
    PUSH_MPS("sing44", 54745, 59708, 8128831.1772, 8.10251e+06, 6, false);
    PUSH_MPS("snp-02-004-104", 126512, 228350, 586803238.6567289, 5.48045e+08, 6, false);
    PUSH_MPS("sorrell3", 169162, 1024, -16.0, -512, 6, false);
    PUSH_MPS("sp150x300d", 450, 600, 69.0, 4.89111, 5, false);
    PUSH_MPS("sp97ar", 1761, 14101, 660705645.759, 6.5256e+08, 6, false);
    PUSH_MPS("sp98ar", 1435, 15085, 529740623.2, 5.24401e+08, 6, false);
    PUSH_MPS("splice1k1", 6505, 3253, -394.0, -1646.78, 6, false);
    PUSH_MPS("square41", 40160, 62234, 15.0, 8.83961, 6, false);
    PUSH_MPS("square47", 61591, 95030, 15.9999999997877, 8.75581, 6, false);
    PUSH_MPS("supportcase10", 165684, 14770, 7.0, 3.38392, 6, false);
    PUSH_MPS("supportcase12", 166781, 799616, -7559.533053817, -7609.53, 6, false);
    PUSH_MPS("supportcase18", 240, 13410, 48.0, 47.1867, 6, false);
    PUSH_MPS("supportcase19", 10713, 1429098, 12677205.99992052, 1.26772e+07, 6, false);
    PUSH_MPS("supportcase26", 870, 436, 1745.123813, 1288.1, 6, false);
    PUSH_MPS("supportcase33", 20489, 20203, -345.0, -552.5, 6, false);
    PUSH_MPS("supportcase40", 38192, 16440, 24256.3122898, 22536.5, 6, false);
    PUSH_MPS("supportcase42", 18439, 19466, 7.75863072227, 7.34356, 6, false);
    PUSH_MPS("supportcase6", 771, 130052, 51906.47737, 45090.7, 6, false);
    PUSH_MPS("supportcase7", 6532, 138844, -1132.223177, -2396.22, 5, false);
    PUSH_MPS("swath1", 884, 6805, 379.07129575, 334.497, 2, false);
    PUSH_MPS("swath3", 884, 6805, 397.76134365, 334.497, 5, false);
    PUSH_MPS("tbfp-network", 2436, 72747, 24.16319444, 23.2615, 5, false);
    PUSH_MPS("thor50dday", 53360, 106261, 40417.0, 4174.49, 6, false);
    PUSH_MPS("timtab1", 171, 397, 764771.99999978, 28694, 6, false);
    PUSH_MPS("tr12-30", 750, 1080, 130595.9999999999, 14210.4, 6, false);
    PUSH_MPS("traininstance2", 15603, 12890, 71820.0, 0, 6, false);
    PUSH_MPS("traininstance6", 12309, 10218, 28290.0, 0, 6, false);
    PUSH_MPS("trento1", 1265, 7687, 5189487.0, 5.18276e+06, 6, false);
    PUSH_MPS("triptim1", 15706, 30055, 22.8680999999999, 22.8681, 3, false);
    PUSH_MPS("uccase12", 121161, 62529, 11507.4050616, 11507.3, 4, false);
    PUSH_MPS("uccase9", 49565, 33242, 10993.131409, 10819.6, 6, false);
    PUSH_MPS("uct-subprob", 1973, 2256, 314.0, 242, 6, false);
    PUSH_MPS("unitcal_7", 48939, 25755, 19635558.244, 1.93876e+07, 5, false);
    PUSH_MPS("var-smallemery-m6j6", 13416, 5608, -149.37501, -156.056, 6, false);
    PUSH_MPS("wachplan", 1553, 3361, -8.0, -9, 5, false);
    // END MIPLIB 2017 TABLE
  } else {
/*
  Load up the problem vector. Note that the row counts here include the
  objective function.
*/
    if (testSwitch<1000001) {
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
    if (testSwitch==1000001) {
      bool anyChosen = false;
      for (int i=21;i<30;i++) {
	if (doThisSet[i])
	  anyChosen = true;
      }
      if (!anyChosen)
	doThisSet[21]=1;
      // read descriptions from file directory/testset
      /* format of file testset
	 comments start with * e.g.
	 * really only name and objective value needed
	 * difficulty should be 21... 29
	 * use -dirmiplib ~/directory++level21-22 to get more difficult ones
	 * if no level then just 21
	 so format of one test case might be
	 value_2631789469.lp 132238 141689 2.6317998e+09 2631789469 21
      */
      FILE * fp = fopen(test1.c_str(),"r");
      if (!fp) {
	printf("Unable to open file %s\n",test1.c_str());
	exit(77);
      }
      // format is name rows columns solution continuous difficulty
      char line[100];
      char name[50];
      int numberRows;
      int numberColumns;
      double objValueI;
      double contValue;
      int difficulty;

      while(fgets(line,100,fp)) {
	if (line[0]=='\n')
	  continue;
	if (line[0]!='*') {
	  if (strchr(line,'\n'))
	    *strchr(line,'\n')='\0';
	  char * where = line;
	  char * space;
	  space = strchr(where,' ');
	  assert(space);
	  *space = '\0';
	  strcpy(name,where);
	  where = space+1;
	  space = strchr(where,' ');
	  assert(space);
	  *space = '\0';
	  numberRows = atoi(where);
	  where = space+1;
	  space = strchr(where,' ');
	  assert(space);
	  *space = '\0';
	  numberColumns = atoi(where);
	  where = space+1;
	  space = strchr(where,' ');
	  assert(space);
	  *space = '\0';
	  objValueI = atof(where);
	  where = space+1;
	  space = strchr(where,' ');
	  if (space) {
	    *space = '\0';
	    contValue = atof(where);
	    where = space+1;
	    difficulty = atoi(where);
	  } else {
	    contValue = atof(where);
	    difficulty = 21;
	  }
	  printf("%s has %d rows, %d columns - objective %g (continuous %g) difficulty %d\n",
		 name,numberRows,numberColumns,objValueI,contValue,difficulty);
	  PUSH_MPS(name,numberRows,numberColumns,objValueI,contValue,difficulty,false);
	} else {
	  printf("%s",line);
	}
      }
      fclose(fp);
    }
  }
#undef PUSH_MPS

  /*
  Normally the problems are executed in order. Define RANDOM_ORDER below to
  randomize.

  #define RANDOM_ORDER
*/
  int which[400];
  int nLoop = static_cast< int >(mpsName.size());
  assert(nLoop <= 400);
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
  // tidy original input queue
  std::deque<std::string> inputQueue = originalInputQueue;
  // for operations after "solve"
  std::deque<std::string> extraQueue;
  // See if "maximize"
  bool maximize = false;
  for (int i = 0; i < inputQueue.size(); i++) {
    if (inputQueue[i] == "-maximize") {
      maximize = true;
    }
  }
  {
    // The queue replayed for each problem is the command line without the
    // directory settings (wherever they are) and without -unitTest itself.
    // Anything after -unitTest is run after each solve. Names are matched
    // as the command line parser matches them: leading dashes and case
    // don't matter. (Leaving -unitTest in would start the test again for
    // every problem.)
    std::deque<std::string> kept;
    size_t unitTestPos = 0;
    for (; unitTestPos < inputQueue.size(); unitTestPos++) {
      const std::string &field = inputQueue[unitTestPos];
      size_t start = field.find_first_not_of('-');
      std::string name = start == std::string::npos ? "" : field.substr(start);
      for (size_t k = 0; k < name.size(); k++)
        name[k] = static_cast<char>(tolower(name[k]));
      if (name == "dirmiplib" || name == "dirsample") {
        unitTestPos++; // and its value
        continue;
      }
      if (name == "unittest")
        break;
      kept.push_back(field);
    }
    if (unitTestPos + 1 < inputQueue.size()) {
      std::cout << "Warning - extra operations after solve" << std::endl;
      // put on extra queue
      for (size_t i = unitTestPos + 1; i < inputQueue.size(); i++) {
        extraQueue.push_back(inputQueue[i]);
        std::cout << inputQueue[i] << " ";
      }
      std::cout << std::endl;
    }
    inputQueue = kept;
  }

//#define CLP_FACTORIZATION_INSTRUMENT
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double timeTakenFac = 0.0;
#endif
  int totalNodes = 0;
  int totalIterations = 0;
  int totalStrong = 0;
  int totalStrongIts = 0;
  int totalStrongFathomed = 0;
  int totalStrongFixed = 0;
  int totalFathomCount = 0;
  int totalExtraNodes = 0;
  int totalExtraIterations = 0;
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
      std::deque<std::string> newInputQueue;
      char replace[100];
      //newArgv[0] = "unitTestCbc";
      if (maximize) {
	solver1.readMps(fn.c_str(), "");
	double objectiveOffset;
	solver1.getDblParam(OsiObjOffset, objectiveOffset);
	solver1.setDblParam(OsiObjOffset, -objectiveOffset);
	int numberColumns = solver1.getNumCols();
	double *array =
	  CoinCopyOfArray(solver1.getObjCoefficients(), numberColumns);
	for (int i = 0; i < numberColumns; i++)
	  array[i] = -array[i];
	solver1.setObjective(array);
	// keep accuracy
	solver1.writeMpsNative("maxVersion.mps",NULL,NULL,2);
	delete[] array;
	fn = "maxVersion.mps";
      }
      newInputQueue.push_back(fn);
      for (int i = 0; i < inputQueue.size(); i++) {
	if (!strstr(inputQueue[i].c_str(),"++")) {
           if (testSwitch >=1000000) {
              // take out dextra3
              if (inputQueue[i] == "dextra3") {
                 continue;
              }
           }
	   newInputQueue.push_back(inputQueue[i]);
	} else {
          //FIXME: This should be changed to use modern C++
          int n = strstr(inputQueue[i].c_str(), "++") - inputQueue[i].c_str();
	  strncpy(replace, inputQueue[i].c_str(), n);
	  const char * mipname = mpsName[m].c_str();
	  int n1 = n;
	  for (int j=0;j<strlen(mipname);j++){
	    replace[n++]=mipname[j];
          }
	  for (int j=n1+2;j<inputQueue[i].length();j++){
             replace[n++]=inputQueue[i].c_str()[j];
          }
	  replace[n] = '\0';
	  newInputQueue.push_back(replace);
	  printf("Replacing %s by %s\n",inputQueue[i].c_str(),replace);
	}
      }
      /*
	Activate the row cut debugger, if requested.
      */
      if (rowCutDebugger[m] == true) {
         newInputQueue.push_back("-debug");
         newInputQueue.push_back("unitTest");
      }
      if (newInputQueue.back()!="-solve")
	newInputQueue.push_back("-solve");
      // add in extra commands
      for (int i = 0; i < extraQueue.size(); i++) {
	if (!strstr(extraQueue[i].c_str(),"++")) {
	   newInputQueue.push_back(extraQueue[i]);
	} else {
          //FIXME: This should be changed to use modern C++
          int n = strstr(extraQueue[i].c_str(), "++") - extraQueue[i].c_str();
	  strncpy(replace, extraQueue[i].c_str(), n);
	  const char * mipname = mpsName[m].c_str();
	  int n1 = n;
	  for (int j=0;j<strlen(mipname);j++){
	    replace[n++]=mipname[j];
          }
	  for (int j=n1+2;j<extraQueue[i].length();j++){
             replace[n++]=extraQueue[i].c_str()[j];
          }
	  replace[n] = '\0';
	  newInputQueue.push_back(replace);
	  printf("Replacing %s by %s\n",extraQueue[i].c_str(),replace);
	}
      }
      model = new CbcModel(solver1);
      CbcMain0(*model, parameters);
      CbcMain1(newInputQueue, *model, parameters, callBack);
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
	      smallest[iRow] = std::min(smallest[iRow], value);
	      largest[iRow] = std::max(largest[iRow], value);
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
	      largest = std::max(largest, above);
	    }
	    if (below < 1.0e12) {
	      largest = std::max(largest, below);
	    }
	    if (rowScale) {
	      double multiplier = rowScale[iRow];
	      above *= multiplier;
	      below *= multiplier;
	    }
	    if (above < 1.0e12) {
	      largestScaled = std::max(largestScaled, above);
	    }
	    if (below < 1.0e12) {
	      largestScaled = std::max(largestScaled, below);
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
	      largest = std::max(largest, above);
	    }
	    if (below < 1.0e12) {
	      largest = std::max(largest, below);
	    }
	    if (columnScale) {
	      double multiplier = 1.0 / columnScale[iColumn];
	      above *= multiplier;
	      below *= multiplier;
	    }
	    if (above < 1.0e12) {
	      largestScaled = std::max(largestScaled, above);
	    }
	    if (below < 1.0e12) {
	      largestScaled = std::max(largestScaled, below);
	    }
	  }
	  std::cout << "Largest (scaled) away from bound " << largestScaled
		    << " unscaled " << largest << std::endl;
#ifdef JJF_ZERO
	  modelC->setDualBound(std::max(1.0001e8,
				       std::min(1000.0 * largestScaled, 1.00001e10)));
#else
	  modelC->setDualBound(std::max(1.0001e9,
				       std::min(1000.0 * largestScaled, 1.0001e10)));
#endif
	}
      } // end clp-specific setup
      /*
	Cut passes: For small models (n < 500) always do 100 passes, if possible
	(-100). For larger models, use minimum drop to stop (100, 20).
      */
      model->setMinimumDrop(std::min(5.0e-2,
				    fabs(model->getMinimizationObjValue()) * 1.0e-3 + 1.0e-4));
      if (std::abs(model->getMaximumCutPassesAtRoot()) <= 100) {
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
    if (maximize)
      objActual = -objActual;
    double objExpect = objValue[m];
    double tolerance = std::min(fabs(objActual), fabs(objExpect));
    tolerance = std::max(1.0e-4, 1.0e-5 * tolerance);
    // The search may stop as soon as it is within the run's own gap
    // tolerances (-allowableGap, -ratioGap; their defaults are 1e-6 and 1e-4,
    // where 2.10's were 0), so a solution worse than the optimum by up to
    // that much is a correct result. (objActual is in the original
    // minimization sense, also with -maximize.)
    double gapTolerance = std::max(model->getAllowableGap(),
      model->getAllowableFractionGap()
        * std::max(fabs(objActual), fabs(objExpect)));
    if (!model->status()) {

      //CoinRelFltEq eq(1.0e-3) ;

      std::cout
        << "cbc_clp (" << mpsName[m] << ") "
        << std::setprecision(10) << objActual;
      if (fabs(objActual - objExpect) < tolerance) {
        std::cout << std::setprecision(dfltPrecision) << "; okay";
        numProbSolved++;
      } else if (objExpect != -1.0e50 && objActual > objExpect
        && objActual - objExpect < tolerance + gapTolerance) {
        std::cout << " (optimum " << objExpect << std::setprecision(dfltPrecision)
                  << "); okay within the gap tolerance";
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
    totalNodes += model->getNodeCount();
    totalIterations += model->getIterationCount();
    const int * strongInfo = model->strongInfo();
    totalStrong += strongInfo[0];
    totalStrongIts += model->numberStrongIterations();
    totalStrongFathomed += strongInfo[2];
    totalStrongFixed += strongInfo[1];
    totalFathomCount += model->getFathomCount();
    totalExtraNodes += model->getExtraNodeCount();
    totalExtraIterations += model->numberExtraIterations();
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
    std::cout
      << " and took " << timeTaken << " seconds." << std::endl;
  } else {
    std::cout
      << " and took " << timeTaken << " seconds." << std::endl;
    std::cout << "cbc_clp Total nodes " << totalNodes << " and " << totalIterations << " iterations - time " << timeTaken << std::endl;
#if 0 // out until statistics fixed
    std::cout << "cbc_clp Total Strong " << totalStrong << " times (" <<
      totalStrongIts << " iterations) - fathomed " <<
      totalStrongFathomed << " - fixed " << totalStrongFixed << std::endl;
    std::cout << "cbc_clp Total Fast node fathoming " << totalFathomCount
	      << " times (" << totalExtraNodes << " nodes, "
	      << totalExtraIterations << " iterations)" << std::endl;
#endif
  }
  if (testSwitch == -2) {
    if (numberFailures || numberOnNodes) {
      std::cout << "****** Unit Test failed." << std::endl;
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
