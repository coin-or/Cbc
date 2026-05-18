// Copyright (C) 2020, COIN-OR Foundation
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file OsiFeatures.hpp
    \brief Defines problem features which can be used to perform algorithm
    parameter recommendation.

    This source contains a list of problem features that can be extracted from a Mixed-Integer Linear
    Program (MIP) from an OsiSolverInterface object. All features are numeric and stored as double.
    All features are extracted in at most O(nz) time, where nz is the number of non-zeros in the
    constraint matrix. Currently 207 features are extracted.

    Some details on algorithm recommendation for CLP can be found in:

    Vilas Boas, M. G.; Santos, H.G.; Merschmann, L.H.C. and Vanden Berghe, G.
    Optimal Decision Trees for the Algorithm Selection Problem:
    Integer Programming Based Approaches. International Transactions in Operational
    Research, DOI 10.1111/itor.12724. 2019.

*/

#include "OsiConfig.h"

class OsiSolverInterface;

/*! List of all features that are extracted */
enum OsiFeature
{
  OFcols = 0,  //< number of columns (variables)
  OFrows,  //< number of rows (constraints)
  OFcolsPerRow, //< cols/rows
  OFequalities, //< number of constraints =
  OFNzEqualities, //< number of non-zeros in equalities
  OFpercEqualities, //< percentage of inequalities
  OFpercNzEqualities, //< percentage of non-zero elements in equalities
  OFinequalities, //< number of constraints <= or >=
  OFNzInequalities, //< number of constraints <= or >=
  OFnz, //< number of non-zero elements in the constraint matrix
  OFdensity, //< density of the constraint matrix, i.e. (nz / (cols*rows))*100

  OFbin, //< number of binary variables
  OFgenInt, //< total number of integer (excluding binaries)
  OFinteger, //< total number of integer (includes binaries) variables
  OFcontinuous, //< number of continuous variables
  OFpercInteger, //< percentage of integer variables
  OFpercBin, //< percentage of binary variables
  OFnUnbounded1, //< number of unbounded variables for one sense (positive or negative)
  OFpercUnbounded1, //< number of unbounded variables for one sense (positive or negative)
  OFnUnbounded2, //< number of unbounded variables for both senses
  OFpercUnbounded2, //< number of unbounded variables for both senses


  /* row types */
  OFrowsPartitioning, //< number of partitioning constraints, e.g. x1 + x2 ... = 1 (binary vars)
  OFpercRowsPartitioning, //< percentage of partitioning constraints, e.g. x1 + x2 ... = 1 (binary vars)
  OFrowsPacking, //< number of packing constraints, e.g. x1 + x2 ... <= 1 (binary vars)
  OFpercRowsPacking, //< percentage of packing constraints, e.g. x1 + x2 ... <= 1 (binary vars)
  OFrowsPartPacking, //< number of partitioning or packing constraints
  OFpercRowsPartPacking, //< percentage of partitioning or packing
  OFrowsCovering, //< number of covering constraints, e.g. x1 + x2 ... >= 1 (binary vars)
  OFpercRowsCovering, //< percentage of covering constraints, e.g. x1 + x2 ... >= 1 (binary vars)
  OFrowsCardinality, //< number of cardinality constraints, e.g. x1 + x2 ... = k, k >= 2 (binary vars)
  OFpercRowsCardinality, //< percentage of cardinality constraints, e.g. x1 + x2 ... = k, k >= 2 (binary vars)
  OFrowsKnapsack, //< number of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, b integer >= 2 (binary vars)
  OFpercRowsKnapsack, //< percentage of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, b integer >= 2 (binary vars)
  OFrowsIntegerKnapsack, //< number of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, cj and b integer, b >= 2 (binary vars)
  OFpercRowsIntegerKnapsack, //< percentage of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, cj and b integer, b >= 2 (binary vars)
  OFrowsInvKnapsack, //< number of invariant knapsack constraints, e.g. x1 + x2 ... <= b, b >= 2 (binary vars)
  OFpercRowsInvKnapsack, //< percentage of invariant knapsack constraints, e.g. x1 + x2 ... <= b, b >= 2 (binary vars)
  OFrowsSingleton, //< number of constraints with only one variable
  OFpercRowsSingleton, //< percentage of constraints with only one variable
  OFrowsAggr, //< number of constraints in the form ax + by = c
  OFpercRowsAggr, //< percentage of constraints in the form ax + by = c
  OFrowsPrec, //< number of constraints in the form ax - ay <= b
  OFpercRowsPrec, //< percentage of constraints in the form ax - ay <= b
  OFrowsVarBnd, //< number of constraints with only one variable
  OFpercRowsVarBnd, //< percentage of constraints with only one variable
  OFrowsBinPacking, //< number of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, cj and b integer, b >= 2, at least one cj >= 2 (binary vars)
  OFpercRowsBinPacking, //< percentage of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, cj and b integer, b >= 2, at least one cj >= 2 (binary vars)
  OFrowsMixedBin, //< constraint that involves binary and continuous variables
  OFpercRowsMixedBin, //< percentage constraint that involves binary and continuous variables
  OFrowsGenInt,  //< constraints with some general integer (not binary) variables
  OFpercRowsGenInt,  //< percentage constraints with some general integer (not binary) variables
  OFrowsFlowBin, //< equality, at least 2 positive and 2 negative coefficients, only bin vars
  OFpercRowsFlowBin, //< percentage equality, at least 2 positive and 2 negative coefficients, only bin vars
  OFrowsFlowMx, //< equality, at least 2 positive and 2 negative coefficients
  OFpercRowsFlowMx, //< equality, at least 2 positive and 2 negative coefficients

  // non zeros distribution per constraint type
  OFnzRowsPartitioning,
  OFnzPercRowsPartitioning,
  OFnzRowsPacking,
  OFnzPercRowsPacking,
  OFnzrowsPartPacking,
  OFnzpercRowsPartPacking,
  OFnzRowsCovering,
  OFnzPercRowsCovering,
  OFnzRowsCardinality,
  OFnzPercRowsCardinality,
  OFnzRowsKnapsack,
  OFnzPercRowsKnapsack,
  OFnzRowsIntegerKnapsack,
  OFnzPercRowsIntegerKnapsack,
  OFnzRowsInvKnapsack,
  OFnzPercRowsInvKnapsack,
  OFnzRowsSingleton,
  OFnzPercRowsSingleton,
  OFnzRowsAggr,
  OFnzPercRowsAggr,
  OFnzRowsPrec,
  OFnzPercRowsPrec,
  OFnzRowsVarBnd,
  OFnzPercRowsVarBnd,
  OFnzRowsBinPacking,
  OFnzPercRowsBinPacking,
  OFnzRowsMixedBin,
  OFnzPercRowsMixedBin,
  OFnzRowsGenInt,
  OFnzPercRowsGenInt,
  OFnzRowsFlowBin,
  OFnzPercRowsFlowBin,
  OFnzRowsFlowMx,
  OFnzPercRowsFlowMx,

  /* statistics constraint matrix */
  OFaMin,
  OFaMax,
  OFaAvg,
  OFaStdDev,
  OFaRatioLSA,
  OFaAllInt,
  OFaPercInt,
  OFaDiffVal,
  OFanShortInts,
  OFapercShortInts,

  /* statistics objective function */
  OFobjMin,
  OFobjMax,
  OFobjAvg,
  OFobjStdDev,
  OFobjRatioLSA,
  OFobjAllInt,
  OFobjPercInt,
  OFobjDiffVal,
  OFobjnShortInts,
  OFobjpercShortInts,

  /* statistics right hand side */
  OFrhsMin,
  OFrhsMax,
  OFrhsAvg,
  OFrhsStdDev,
  OFrhsRatioLSA,
  OFrhsAllInt,
  OFrhsPercInt,
  OFrhsDiffVal,
  OFrhsnShortInts,
  OFrhspercShortInts,

  /* statistics non-zero distribution rows */
  OFrowNzMin,
  OFrowNzMax,
  OFrowNzAvg,
  OFrowNzStdDev,

  /* statistics non-zero distribution cols */
  OFcolNzMin,
  OFcolNzMax,
  OFcolNzAvg,
  OFcolNzStdDev,

  // constraints with nz less or equal
  OFrowsLess4Nz,
  OFrowsLess8Nz,
  OFrowsLess16Nz,
  OFrowsLess32Nz,
  OFrowsLess64Nz,
  OFrowsLess128Nz,
  OFrowsLess256Nz,
  OFrowsLess512Nz,
  OFrowsLess1024Nz,
  OFpercRowsLess4Nz,
  OFpercRowsLess8Nz,
  OFpercRowsLess16Nz,
  OFpercRowsLess32Nz,
  OFpercRowsLess64Nz,
  OFpercRowsLess128Nz,
  OFpercRowsLess256Nz,
  OFpercRowsLess512Nz,
  OFpercRowsLess1024Nz,

  // constraints nz at least
  OFrowsLeast4Nz,
  OFrowsLeast8Nz,
  OFrowsLeast16Nz,
  OFrowsLeast32Nz,
  OFrowsLeast64Nz,
  OFrowsLeast128Nz,
  OFrowsLeast256Nz,
  OFrowsLeast512Nz,
  OFrowsLeast1024Nz,
  OFrowsLeast2048Nz,
  OFrowsLeast4096Nz,
  OFpercRowsLeast4Nz,
  OFpercRowsLeast8Nz,
  OFpercRowsLeast16Nz,
  OFpercRowsLeast32Nz,
  OFpercRowsLeast64Nz,
  OFpercRowsLeast128Nz,
  OFpercRowsLeast256Nz,
  OFpercRowsLeast512Nz,
  OFpercRowsLeast1024Nz,
  OFpercRowsLeast2048Nz,
  OFpercRowsLeast4096Nz,

  // constraints with nz less or equal
  OFcolsLess4Nz,
  OFcolsLess8Nz,
  OFcolsLess16Nz,
  OFcolsLess32Nz,
  OFcolsLess64Nz,
  OFcolsLess128Nz,
  OFcolsLess256Nz,
  OFcolsLess512Nz,
  OFcolsLess1024Nz,
  OFpercColsLess4Nz,
  OFpercColsLess8Nz,
  OFpercColsLess16Nz,
  OFpercColsLess32Nz,
  OFpercColsLess64Nz,
  OFpercColsLess128Nz,
  OFpercColsLess256Nz,
  OFpercColsLess512Nz,
  OFpercColsLess1024Nz,

  // constraints nz at least
  OFcolsLeast4Nz,
  OFcolsLeast8Nz,
  OFcolsLeast16Nz,
  OFcolsLeast32Nz,
  OFcolsLeast64Nz,
  OFcolsLeast128Nz,
  OFcolsLeast256Nz,
  OFcolsLeast512Nz,
  OFcolsLeast1024Nz,
  OFcolsLeast2048Nz,
  OFcolsLeast4096Nz,

  OFpercColsLeast4Nz,
  OFpercColsLeast8Nz,
  OFpercColsLeast16Nz,
  OFpercColsLeast32Nz,
  OFpercColsLeast64Nz,
  OFpercColsLeast128Nz,
  OFpercColsLeast256Nz,
  OFpercColsLeast512Nz,
  OFpercColsLeast1024Nz,
  OFpercColsLeast2048Nz,
  OFpercColsLeast4096Nz,

  OFCount
};

class OSILIB_EXPORT OsiFeatures {
public:
  /** @brief number of features */
  static int n;

  /** @brief name of the i-th feature */
  static const char *name(int i);

  /** @brief name of an specific feature */
  static const char *name( const OsiFeature of );

  /** @brief computes all feature values, the size of this vector should be at least OFCount */
  static void compute(double *features, OsiSolverInterface *solver);
};


/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
