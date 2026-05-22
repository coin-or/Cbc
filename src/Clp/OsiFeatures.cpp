#include "OsiSolverInterface.hpp"
#include <limits>
#include <algorithm>
#include <cmath>
#include <vector>
#include "OsiFeatures.hpp"

int OsiFeatures::n = OFCount;

bool dbl_equal( const double v1, const double v2 ) {
    return fabs(v1-v2) <= 1e-16;
}

bool intVal( const double val ) {
    return fabs(val - round(val)) <= 1e-16;
}

// max size for features
#define STR_SIZE 64

#if __cplusplus <= 199711L  // no C++ 11
#include <set>
typedef std::set<double> uSetD;
#else
#include <unordered_set>
typedef std::unordered_set<double> uSetD;
#endif

// store at most this number of different values
#define MAX_DIFF_VALUES 4096

class Summary {
    public:
        Summary() :
            minV( std::numeric_limits<double>::max() ),
            maxV( std::numeric_limits<double>::lowest() ),
            minAbsV( std::numeric_limits<double>::max() ),
            maxAbsV( std::numeric_limits<double>::lowest() ),
            avg(0.0),
            ratioLSA(0.0),
            percIntEl(0.0),
            summV( 0.0 ),
            nEl(0),
            intEl(0),
            allIntEl(0),
            nPosVal(0),
            nNegVal(0),
            nShortInt(0),
            lastV(COIN_DBL_MAX),
            nDiffValues(0)
    { }

        void add( double val ) {
            this->minV = std::min(val, minV);
            this->maxV = std::max(val, maxV);
            double absv = std::abs(val);
            if (absv >= 1e-16) {
                this->minAbsV = std::min(absv, minAbsV);
                this->maxAbsV = std::max(absv, maxAbsV);
            }
            this->summV += val;
            if (intVal(val)) {
                this->intEl++;
                if (absv < 65535)
                    nShortInt++;
            }
            if (val >= 1e-16)
                this->nPosVal++;
            else {
                if (val <= -1e-16)
                    this->nNegVal++;
            }

            if (val != lastV) {
                lastV = val;
                if (values.size()<MAX_DIFF_VALUES)
                    values.insert(val);
            }

            this->nEl++;
        }

        void finish() {
            this->avg = this->summV / ((long double) this->nEl);
            if (this->minAbsV != 0.0)
                this->ratioLSA = this->maxAbsV / this->minAbsV;

            if (this->intEl == this->nEl)
                this->allIntEl = 1;
            else
                this->allIntEl = 0;

            this->nDiffValues = this->values.size();

            this->percIntEl = (((double)this->intEl) / ((double) this->nEl)) * 100.0;
        }

        double minV;
        double maxV;
        double minAbsV;
        double maxAbsV;
        double avg;
        double ratioLSA;
        double percIntEl;
        long double summV;
        unsigned int nEl;
        unsigned int intEl;
        unsigned int allIntEl;
        unsigned int nPosVal;
        unsigned int nNegVal;
        unsigned int nShortInt;
        uSetD values;
        double lastV;
        unsigned int nDiffValues;
};

const static char feat_names[OFCount][STR_SIZE] = {
    "cols",
    "rows",
    "colsPerRow",
    "equalities",
    "nzEqualities",
    "percEqualities",
    "percNzEqualities",
    "inequalities",
    "nzInequalities",
    "nz",
    "density",

    "bin",
    "genInt",
    "integer",
    "continuous",
    "percInteger",
    "percBin",
    "nUnbounded1",
    "percUnbounded1",
    "nUnbounded2",
    "percUnbounded2",

    "rPartitioning",
    "rPercPartitioning",
    "rPacking",
    "rPercPacking",
    "rPartPacking",
    "rPercRowsPartPacking",
    "rCovering",
    "rPercCovering",
    "rCardinality",
    "rPercCardinality",
    "rKnapsack",
    "rPercKnapsack",
    "rIntegerKnapsack",
    "rPercIntegerKnapsack",
    "rInvKnapsack",
    "rPercInvKnapsack",
    "rSingleton",
    "rPercSingleton",
    "rAggre",
    "rPercAggre",
    "rPrec",
    "rPercPrec",
    "rVarBnd",
    "rPercVarBnd",
    "rBinPacking",
    "rPercBinPacking",
    "rMixedBin",
    "rPercMixedBin",
    "rGenInt",
    "rPercGenInt",
    "rFlowBin",
    "rPercFlowBin",
    "rFlowMx",
    "rPercFlowMx",

    "rNzRowsPartitioning",
    "rNzPercRowsPartitioning",
    "rNzRowsPacking",
    "rNzPercRowsPacking",
    "rNzrowsPartPacking",
    "rNzpercRowsPartPacking",
    "rNzRowsCovering",
    "rNzPercRowsCovering",
    "rNzRowsCardinality",
    "rNzPercRowsCardinality",
    "rNzRowsKnapsack",
    "rNzPercRowsKnapsack",
    "rNzRowsIntegerKnapsack",
    "rNzPercRowsIntegerKnapsack",
    "rNzRowsInvKnapsack",
    "rNzPercRowsInvKnapsack",
    "rNzRowsSingleton",
    "rNzPercRowsSingleton",
    "rNzRowsAggr",
    "rNzPercRowsAggr",
    "rNzRowsPrec",
    "rNzPercRowsPrec",
    "rNzRowsVarBnd",
    "rNzPercRowsVarBnd",
    "rNzRowsBinPacking",
    "rNzPercRowsBinPacking",
    "rNzRowsMixedBin",
    "rNzPercRowsMixedBin",
    "rNzRowsGenInt",
    "rNzPercRowsGenInt",
    "rNzRowsFlowBin",
    "rNzPercRowsFlowBin",
    "rNzRowsFlowMx",
    "rNzPercRowsFlowMx",

    "aMin",
    "aMax",
    "aAvg",
    "aStdDev",
    "aRatioLSA",
    "aAllInt",
    "aPercInt",
    "aDiffVal",
    "anShortInts",
    "apercShortInts",

    "objMin",
    "objMax",
    "objAvg",
    "objStdDev",
    "objRatioLSA",
    "objAllInt",
    "objPercInt",
    "objDiffVal",
    "objnShortInts",
    "objpercShortInts",

    "rhsMin",
    "rhsMax",
    "rhsAvg",
    "rhsStdDev",
    "rhsRatioLSA",
    "rhsAllInt",
    "rhsPercInt",
    "rhsDiffVal",
    "rhsnShortInts",
    "rhspercShortInts",

    "rowNzMin",
    "rowNzMax",
    "rowNzAvg",
    "rowNzStdDev",

    "colNzMin",
    "colNzMax",
    "colNzAvg",
    "colNzStdDev",

    "rowsLess4Nz",
    "rowsLess8Nz",
    "rowsLess16Nz",
    "rowsLess32Nz",
    "rowsLess64Nz",
    "rowsLess128Nz",
    "rowsLess256Nz",
    "rowsLess512Nz",
    "rowsLess1024Nz",
    "percRowsLess4Nz",
    "percRowsLess8Nz",
    "percRowsLess16Nz",
    "percRowsLess32Nz",
    "percRowsLess64Nz",
    "percRowsLess128Nz",
    "percRowsLess256Nz",
    "percRowsLess512Nz",
    "percRowsLess1024Nz",
    "rowsLeast4Nz",
    "rowsLeast8Nz",
    "rowsLeast16Nz",
    "rowsLeast32Nz",
    "rowsLeast64Nz",
    "rowsLeast128Nz",
    "rowsLeast256Nz",
    "rowsLeast512Nz",
    "rowsLeast1024Nz",
    "rowsLeast2048Nz",
    "rowsLeast4096Nz",
    "percRowsLeast4Nz",
    "percRowsLeast8Nz",
    "percRowsLeast16Nz",
    "percRowsLeast32Nz",
    "percRowsLeast64Nz",
    "percRowsLeast128Nz",
    "percRowsLeast256Nz",
    "percRowsLeast512Nz",
    "percRowsLeast1024Nz",
    "percRowsLeast2048Nz",
    "percRowsLeast4096Nz",

    "colsLess4Nz",
    "colsLess8Nz",
    "colsLess16Nz",
    "colsLess32Nz",
    "colsLess64Nz",
    "colsLess128Nz",
    "colsLess256Nz",
    "colsLess512Nz",
    "colsLess1024Nz",
    "percColsLess4Nz",
    "percColsLess8Nz",
    "percColsLess16Nz",
    "percColsLess32Nz",
    "percColsLess64Nz",
    "percColsLess128Nz",
    "percColsLess256Nz",
    "percColsLess512Nz",
    "percColsLess1024Nz",

    "colsLeast4Nz",
    "colsLeast8Nz",
    "colsLeast16Nz",
    "colsLeast32Nz",
    "colsLeast64Nz",
    "colsLeast128Nz",
    "colsLeast256Nz",
    "colsLeast512Nz",
    "colsLeast1024Nz",
    "colsLeast2048Nz",
    "colsLeast4096Nz",

    "percColsLeast4Nz",
    "percColsLeast8Nz",
    "perccolsLeast16Nz",
    "perccolsLeast32Nz",
    "perccolsLeast64Nz",
    "perccolsLeast128Nz",
    "perccolsLeast256Nz",
    "perccolsLeast512Nz",
    "perccolsLeast1024Nz",
    "perccolsLeast2048Nz",
    "perccolsLeast4096Nz",
};

const char *OsiFeatures::name( int i ) {
    return feat_names[i];
}

const char *OsiFeatures::name( const OsiFeature of ) {
    return OsiFeatures::name((int) of);
}

double std_dev( const int *el, const double avg, int n ) {
    long double sum = 0.0;

    for ( int i=0 ; (i<n) ; ++i ) {
        const double diff = (((double)el[i])-avg);
        sum += diff * diff;
    }
    sum /= ((long double)n);

    return (double) sqrtl(sum);
}



double std_dev( const double *el, const double avg, int n ) {
    long double sum = 0.0;

    for ( int i=0 ; (i<n) ; ++i ) {
        const double diff = (el[i]-avg);
        sum += diff * diff;
    }
    sum /= ((long double)n);

    return (double) sqrtl(sum);
}

void OsiFeatures::compute(double *features, OsiSolverInterface *solver) {
    // initializing
    for ( int i=0 ; i<OsiFeatures::n ; ++i  )
        features[i] = 0.0;

    features[OFcols] = solver->getNumCols();
    features[OFrows] = solver->getNumRows();
    features[OFcolsPerRow] = ((double)(solver->getNumCols())) / ((double)(solver->getNumRows()));
    features[OFnz] = solver->getNumElements();
    features[OFdensity] = (((long double)solver->getNumElements()) / (((long double)solver->getNumCols())*((long double)solver->getNumRows()))) * ((long double) 100.0);

    // Precompute per-column type to avoid O(nz) virtual calls inside the row loop.
    // colType[j]: 0=continuous, 1=binary, 2=general integer
    const int nCols = solver->getNumCols();
    const double *colLBpre = solver->getColLower();
    const double *colUBpre = solver->getColUpper();
    std::vector<signed char> colType(nCols, 0);
    for (int j = 0; j < nCols; ++j) {
        if (solver->isInteger(j)) {
            const double lb = colLBpre[j], ub = colUBpre[j];
            if ((ub == 1.0 || ub == 0.0) && (lb == 0.0 || lb == 1.0))
                colType[j] = 1; // binary
            else
                colType[j] = 2; // general integer
        }
    }

    Summary aSumm, rhsSumm, objSumm, rowNzSumm, colNzSumm;

    /* going though all rows */
    unsigned int nRows = solver->getNumRows();
    const double *rhs = solver->getRightHandSide();
    const char *sense = solver->getRowSense();
    const CoinPackedMatrix *cpmRow =  solver->getMatrixByRow();
    for ( unsigned int row=0 ; (row<nRows) ; ++row ) {
        const int nzRow = cpmRow->getVectorLengths()[row];
        const CoinBigIndex *starts = cpmRow->getVectorStarts();
        const int *ridx = cpmRow->getIndices() + starts[row];
        const double *rcoef = cpmRow->getElements() + starts[row];

        rowNzSumm.add(nzRow);

        Summary summRow;

        int nBinRow = 0, nContRow = 0, nIntRow = 0;
        for (int j = 0 ; (j < nzRow) ; ++j) {
            aSumm.add(rcoef[j]);
            summRow.add(rcoef[j]);

            switch (colType[ridx[j]]) {
                case 1: nBinRow++;  break;
                case 2: nIntRow++;  break;
                default: nContRow++; break;
            }
        }

        bool rowBin = (nBinRow == nzRow);
        bool intCoefs = (summRow.intEl == summRow.nEl);

        switch (nzRow) {
            case 1:
                features[OFrowsSingleton]++;
                features[OFnzRowsSingleton]++;
                break;
            case 2:
                if (sense[row] == 'E') {
                    features[OFrowsAggr]++;
                    features[OFnzRowsAggr] += nzRow;
                }

                if (nBinRow == 1) {
                    features[OFrowsVarBnd]++;
                    features[OFnzRowsVarBnd] += nzRow;
                }

                if ( nBinRow%2 == 0 && nContRow%2 == 0) // vars of the same type
                    if ( summRow.nNegVal == 1 && summRow.nPosVal == 1
                            && dbl_equal(summRow.minV, summRow.maxV) ) {
                        features[OFrowsPrec]++;
                        features[OFnzRowsPrec] += nzRow;
                    }
                break;
        }

        // constraint types with only binary variables
        const double minaRow = summRow.minV;
        const double maxaRow = summRow.maxV;

        if (rowBin) {
            // pack, part and cov
            // Handle both canonical form (all +1 coefs) and negated form (all -1 coefs,
            // flipped sense/RHS) — e.g. sorrell3 stores x_i+x_j<=1 as -x_i-x_j>=-1.
            const bool allPlusOnes = dbl_equal(minaRow, 1.0) && dbl_equal(maxaRow, 1.0);
            const bool allMinusOnes = dbl_equal(minaRow, -1.0) && dbl_equal(maxaRow, -1.0);
            // Normalise negated rows: flip sense and negate rhs so the rest of the logic
            // only needs to handle the canonical (+1 coefs) form.
            const char effectiveSense = allMinusOnes
                ? (sense[row] == 'L' ? 'G' : (sense[row] == 'G' ? 'L' : sense[row]))
                : sense[row];
            const double effectiveRhs = allMinusOnes ? -rhs[row] : rhs[row];

            if ( allPlusOnes || allMinusOnes ) {
                if (dbl_equal(effectiveRhs, 1.0)) {
                    switch (effectiveSense) {
                        case 'E':
                            features[OFrowsPartitioning]++;
                            features[OFnzRowsPartitioning] += nzRow;
                            features[OFnzrowsPartPacking] += nzRow;
                            break;
                        case 'G':
                            features[OFrowsCovering]++;
                            features[OFnzRowsCovering] += nzRow;
                            break;
                        case 'L':
                            features[OFrowsPacking]++;
                            features[OFnzRowsPacking] += nzRow;
                            features[OFnzrowsPartPacking] += nzRow;
                            break;
                    }
                } // rhs 1.0
                else {
                    if (effectiveRhs >= 1.99) {
                        switch (effectiveSense) {
                            case 'E':
                                features[OFrowsCardinality]++;
                                features[OFnzRowsCardinality] += nzRow;
                                break;
                            case 'L':
                                features[OFrowsInvKnapsack]++;
                                features[OFnzRowsInvKnapsack] += nzRow;
                                break;
                        }
                    } // rhs >= 2
                } // rhs != 1
            } // only ones (or minus-ones) LHS as coefs
            else {
                if (rhs[row] >= 1.1 ) {
                    if ( (maxaRow - minaRow >= 0.1) && (summRow.nNegVal == 0) ) { // different weights
                        features[OFrowsKnapsack]++;
                        features[OFnzRowsKnapsack] += nzRow;
                        if (intCoefs) {
                            features[OFrowsIntegerKnapsack]++;
                            features[OFnzRowsIntegerKnapsack] += nzRow;
                        }
                    }
                    if (summRow.nNegVal == 1 && nzRow >= 2) {
                        features[OFrowsBinPacking]++;
                        features[OFnzRowsBinPacking] += nzRow;
                    }
                }
            }

            if (summRow.nNegVal >= 2 && summRow.nPosVal >= 2 && sense[row]=='E') {
                features[OFrowsFlowBin]++;
                features[OFnzRowsFlowBin] += nzRow;
            }
        } // only binary vars
        else
        {
            if (summRow.nNegVal >= 2 && summRow.nPosVal >= 2 && sense[row]=='E') {
                features[OFrowsFlowMx]++;
                features[OFnzRowsFlowMx] += nzRow;
            }
            if (nContRow>0) {
                features[OFrowsMixedBin]++;
                features[OFnzRowsMixedBin] += nzRow;
            }
            if (nIntRow) {
                features[OFrowsGenInt]++;
                features[OFnzRowsGenInt] += nzRow;
            }
        }

        switch (sense[row]) {
            case 'E':
                features[OFequalities]++;
                features[OFNzEqualities] += nzRow;

                break;
            default:
                // inequalities
                features[OFinequalities]++;
                features[OFNzInequalities] += nzRow;
                break;
        }

        rhsSumm.add( rhs[row] );

        if (nzRow <= 1024) {
            features[OFrowsLess1024Nz]++;
            if (nzRow <= 512) {
                features[OFrowsLess512Nz]++;
                if (nzRow <= 256) {
                    features[OFrowsLess256Nz]++;
                    if (nzRow <= 128) {
                        features[OFrowsLess128Nz]++;
                        if (nzRow <= 64) {
                            features[OFrowsLess64Nz]++;
                            if (nzRow <= 32) {
                                features[OFrowsLess32Nz]++;
                                if (nzRow <= 16) {
                                    features[OFrowsLess16Nz]++;
                                    if (nzRow <= 8) {
                                        features[OFrowsLess8Nz]++;
                                        if (nzRow <= 4) {
                                            features[OFrowsLess4Nz]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } // checking small rows

        if (nzRow >= 4) {
            features[OFrowsLeast4Nz]++;
            if (nzRow >= 8) {
                features[OFrowsLeast8Nz]++;
                if (nzRow >= 16) {
                    features[OFrowsLeast16Nz]++;
                    if (nzRow >= 32) {
                        features[OFrowsLeast32Nz]++;
                        if (nzRow >= 64) {
                            features[OFrowsLeast64Nz]++;
                            if (nzRow >= 128) {
                                features[OFrowsLeast128Nz]++;
                                if (nzRow >= 256) {
                                    features[OFrowsLeast256Nz]++;
                                    if (nzRow >= 512) {
                                        features[OFrowsLeast512Nz]++;
                                        if (nzRow >= 1024) {
                                            features[OFrowsLeast1024Nz]++;
                                            if (nzRow >= 2048) {
                                                features[OFrowsLeast2048Nz]++;
                                                if (nzRow >= 4096) {
                                                    features[OFrowsLeast4096Nz]++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } // larger rows
    } // rows

    features[OFpercRowsPartitioning] = (features[OFrowsPartitioning] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsPacking] = (features[OFrowsPacking] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsPartPacking] = (features[OFrowsPartPacking] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsCovering] = (features[OFrowsCovering] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsCardinality] = (features[OFrowsCardinality] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsIntegerKnapsack] = (features[OFrowsIntegerKnapsack] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsInvKnapsack] = (features[OFrowsInvKnapsack] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsSingleton] = (features[OFrowsSingleton] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsAggr] = (features[OFrowsAggr] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsPrec] = (features[OFrowsPrec] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsVarBnd] = (features[OFrowsVarBnd] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsBinPacking] = (features[OFrowsBinPacking] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsMixedBin] = (features[OFrowsMixedBin] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsGenInt] = (features[OFrowsGenInt] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsFlowBin] = (features[OFpercRowsFlowBin] / (double)solver->getNumRows())*100.0;
    features[OFpercRowsFlowMx] = (features[OFpercRowsFlowMx] / (double)solver->getNumRows())*100.0;

    features[OFnzPercRowsPartitioning] = (features[OFnzRowsPartitioning] / features[OFnz])*100.0;
    features[OFnzpercRowsPartPacking] = (features[OFnzrowsPartPacking] / features[OFnz])*100.0;
    features[OFnzPercRowsPacking] = (features[OFnzRowsPacking] / features[OFnz])*100.0;
    features[OFnzPercRowsCovering] = (features[OFnzRowsCovering] / features[OFnz])*100.0;
    features[OFnzPercRowsCardinality] = (features[OFnzRowsCardinality] / features[OFnz])*100.0;
    features[OFnzPercRowsKnapsack] = (features[OFnzRowsKnapsack] / features[OFnz])*100.0;
    features[OFnzPercRowsIntegerKnapsack] = (features[OFnzRowsIntegerKnapsack] / features[OFnz])*100.0;
    features[OFnzPercRowsInvKnapsack] = (features[OFnzRowsInvKnapsack] / features[OFnz])*100.0;
    features[OFnzPercRowsSingleton] = (features[OFnzRowsSingleton] / features[OFnz])*100.0;
    features[OFnzPercRowsAggr] = (features[OFnzRowsAggr] / features[OFnz])*100.0;
    features[OFnzPercRowsPrec] = (features[OFnzRowsPrec] / features[OFnz])*100.0;
    features[OFnzPercRowsVarBnd] = (features[OFnzRowsVarBnd] / features[OFnz])*100.0;
    features[OFnzPercRowsBinPacking] = (features[OFnzRowsBinPacking] / features[OFnz])*100.0;
    features[OFnzPercRowsMixedBin] = (features[OFnzRowsMixedBin] / features[OFnz])*100.0;
    features[OFnzPercRowsGenInt] = (features[OFnzRowsGenInt] / features[OFnz])*100.0;
    features[OFnzPercRowsFlowBin] = (features[OFnzRowsFlowBin] / features[OFnz])*100.0;
    features[OFnzPercRowsFlowMx] = (features[OFnzRowsFlowMx] / features[OFnz])*100.0;

    features[OFpercNzEqualities] = (features[OFNzEqualities] / features[OFnz])*100.0;

    aSumm.finish();
    rhsSumm.finish();
    rowNzSumm.finish();

    // cols
    const CoinPackedMatrix *cpmCol =  solver->getMatrixByCol();
    const double *obj = solver->getObjCoefficients();
    const double *colLB = solver->getColLower();
    const double *colUB = solver->getColUpper();
    for ( int j=0 ; (j<solver->getNumCols()) ; ++j ) {
        objSumm.add(obj[j]);

        switch (colType[j]) {
            case 1:
                features[OFbin]++;
                features[OFinteger]++;
                break;
            case 2:
                features[OFgenInt]++;
                features[OFinteger]++;
                break;
            default:
                features[OFcontinuous]++;
                break;
        }

        if (colUB[j] == COIN_DBL_MAX && colLB[j] == -COIN_DBL_MAX) {
            features[OFnUnbounded2]++;
        } else {
            if (colUB[j] == COIN_DBL_MAX || colLB[j] == -COIN_DBL_MAX) {
                features[OFnUnbounded1]++;
            }
        }

        const int nzCol = cpmCol->getVectorLengths()[j];

        colNzSumm.add(nzCol);
        //const CoinBigIndex *starts = cpmCol->getVectorStarts();
        //const int *ridx = cpmCol->getIndices() + starts[j];
        //const double *rcoef = cpmCol->getElements() + starts[j];

        if (nzCol <= 1024) {
            features[OFcolsLess1024Nz]++;
            if (nzCol <= 512) {
                features[OFcolsLess512Nz]++;
                if (nzCol <= 256) {
                    features[OFcolsLess256Nz]++;
                    if (nzCol <= 128) {
                        features[OFcolsLess128Nz]++;
                        if (nzCol <= 64) {
                            features[OFcolsLess64Nz]++;
                            if (nzCol <= 32) {
                                features[OFcolsLess32Nz]++;
                                if (nzCol <= 16) {
                                    features[OFcolsLess16Nz]++;
                                    if (nzCol <= 8) {
                                        features[OFcolsLess8Nz]++;
                                        if (nzCol <= 4) {
                                            features[OFcolsLess4Nz]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } // checking small cols

        if (nzCol >= 4) {
            features[OFcolsLeast4Nz]++;
            if (nzCol >= 8) {
                features[OFcolsLeast8Nz]++;
                if (nzCol >= 16) {
                    features[OFcolsLeast16Nz]++;
                    if (nzCol >= 32) {
                        features[OFcolsLeast32Nz]++;
                        if (nzCol >= 64) {
                            features[OFcolsLeast64Nz]++;
                            if (nzCol >= 128) {
                                features[OFcolsLeast128Nz]++;
                                if (nzCol >= 256) {
                                    features[OFcolsLeast256Nz]++;
                                    if (nzCol >= 512) {
                                        features[OFcolsLeast512Nz]++;
                                        if (nzCol >= 1024) {
                                            features[OFcolsLeast1024Nz]++;
                                            if (nzCol >= 2048) {
                                                features[OFcolsLeast2048Nz]++;
                                                if (nzCol >= 4096) {
                                                    features[OFcolsLeast4096Nz]++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } // larger cols

    } // all variables
    objSumm.finish();
    features[OFpercInteger] = (features[OFinteger] / ((double) solver->getNumCols())) * 100.0;
    features[OFpercBin] = (features[OFbin] / ((double) solver->getNumCols())) * 100.0;
    features[OFpercUnbounded1] = (features[OFnUnbounded1] / ((double) solver->getNumCols())) * 100.0;
    features[OFpercUnbounded2] = (features[OFnUnbounded2] / ((double) solver->getNumCols())) * 100.0;

    features[OFaMin] = aSumm.minV;
    features[OFaMax] = aSumm.maxV;
    features[OFaAvg] = aSumm.avg;
    features[OFaStdDev] = std_dev( cpmRow->getElements(), aSumm.avg, solver->getNumElements() );
    features[OFaRatioLSA] = aSumm.ratioLSA;
    features[OFaAllInt] = aSumm.allIntEl;
    features[OFaPercInt] = aSumm.percIntEl;
    features[OFaDiffVal] = aSumm.nDiffValues;
    features[OFanShortInts] = aSumm.nShortInt;
    features[OFapercShortInts] = (features[OFanShortInts] / ((double) solver->getNumElements()))*100.0;

    features[OFobjMin] = objSumm.minV;
    features[OFobjMax] = objSumm.maxV;
    features[OFobjAvg] = objSumm.avg;
    features[OFobjStdDev] = std_dev( solver->getObjCoefficients(), objSumm.avg, solver->getNumCols() );
    features[OFobjRatioLSA] = objSumm.ratioLSA;
    features[OFobjAllInt] = objSumm.allIntEl;
    features[OFobjPercInt] = objSumm.percIntEl;
    features[OFobjDiffVal] = objSumm.nDiffValues;
    features[OFobjnShortInts] = objSumm.nShortInt;
    features[OFobjpercShortInts] = (features[OFobjnShortInts] / ((double)objSumm.nEl))*100.0;

    features[OFrhsMin] = rhsSumm.minV;
    features[OFrhsMax] = rhsSumm.maxV;
    features[OFrhsAvg] = rhsSumm.avg;
    features[OFrhsStdDev] = std_dev( solver->getRightHandSide(), rhsSumm.avg, solver->getNumRows() );
    features[OFrhsRatioLSA] = rhsSumm.ratioLSA;
    features[OFrhsAllInt] = rhsSumm.allIntEl;
    features[OFrhsPercInt] = rhsSumm.percIntEl;
    features[OFrhsDiffVal] = rhsSumm.nDiffValues;
    features[OFrhsnShortInts] = rhsSumm.nShortInt;
    features[OFrhspercShortInts] = (features[OFrhsnShortInts] / ((double)rhsSumm.nEl))*100.0;

    features[OFrowNzMin] = rowNzSumm.minV;
    features[OFrowNzMax] = rowNzSumm.maxV;
    features[OFrowNzAvg] = rowNzSumm.avg;
    features[OFrowNzStdDev] = std_dev(cpmRow->getVectorLengths(), rowNzSumm.avg, solver->getNumRows() );

    features[OFcolNzMin] = colNzSumm.minV;
    features[OFcolNzMax] = colNzSumm.maxV;
    features[OFcolNzAvg] = colNzSumm.avg;
    features[OFcolNzStdDev] = std_dev(cpmCol->getVectorLengths(), colNzSumm.avg, solver->getNumCols() );

    double dnRows = solver->getNumRows();
    features[OFpercRowsLess4Nz] = (features[OFrowsLess4Nz] / dnRows)*100.0;
    features[OFpercRowsLess8Nz] = (features[OFrowsLess8Nz] / dnRows)*100.0;
    features[OFpercRowsLess16Nz] = (features[OFrowsLess16Nz] / dnRows)*100.0;
    features[OFpercRowsLess32Nz] = (features[OFrowsLess32Nz] / dnRows)*100.0;
    features[OFpercRowsLess64Nz] = (features[OFrowsLess64Nz] / dnRows)*100.0;
    features[OFpercRowsLess128Nz] = (features[OFrowsLess128Nz] / dnRows)*100.0;
    features[OFpercRowsLess256Nz] = (features[OFrowsLess256Nz] / dnRows)*100.0;
    features[OFpercRowsLess512Nz] = (features[OFrowsLess512Nz] / dnRows)*100.0;
    features[OFpercRowsLess1024Nz] = (features[OFrowsLess1024Nz] / dnRows)*100.0;

    features[OFpercRowsLeast4Nz] = (features[OFrowsLeast4Nz] / dnRows)*100.0;
    features[OFpercRowsLeast8Nz] = (features[OFrowsLeast8Nz] / dnRows)*100.0;
    features[OFpercRowsLeast16Nz] = (features[OFrowsLeast16Nz] / dnRows)*100.0;
    features[OFpercRowsLeast32Nz] = (features[OFrowsLeast32Nz] / dnRows)*100.0;
    features[OFpercRowsLeast64Nz] = (features[OFrowsLeast64Nz] / dnRows)*100.0;
    features[OFpercRowsLeast128Nz] = (features[OFrowsLeast128Nz] / dnRows)*100.0;
    features[OFpercRowsLeast256Nz] = (features[OFrowsLeast256Nz] / dnRows)*100.0;
    features[OFpercRowsLeast512Nz] = (features[OFrowsLeast512Nz] / dnRows)*100.0;
    features[OFpercRowsLeast1024Nz] = (features[OFrowsLeast1024Nz] / dnRows)*100.0;
    features[OFpercRowsLeast2048Nz] = (features[OFrowsLeast2048Nz] / dnRows)*100.0;
    features[OFpercRowsLeast4096Nz] = (features[OFrowsLeast4096Nz] / dnRows)*100.0;

    double dnCols = solver->getNumCols();
    features[OFpercColsLess4Nz] = (features[OFcolsLess4Nz] / dnCols)*100.0;
    features[OFpercColsLess8Nz] = (features[OFcolsLess8Nz] / dnCols)*100.0;
    features[OFpercColsLess16Nz] = (features[OFcolsLess16Nz] / dnCols)*100.0;
    features[OFpercColsLess32Nz] = (features[OFcolsLess32Nz] / dnCols)*100.0;
    features[OFpercColsLess64Nz] = (features[OFcolsLess64Nz] / dnCols)*100.0;
    features[OFpercColsLess128Nz] = (features[OFcolsLess128Nz] / dnCols)*100.0;
    features[OFpercColsLess256Nz] = (features[OFcolsLess256Nz] / dnCols)*100.0;
    features[OFpercColsLess512Nz] = (features[OFcolsLess512Nz] / dnCols)*100.0;
    features[OFpercColsLess1024Nz] = (features[OFcolsLess1024Nz] / dnCols)*100.0;

    features[OFpercColsLeast4Nz] = (features[OFcolsLeast4Nz] / dnCols)*100.0;
    features[OFpercColsLeast8Nz] = (features[OFcolsLeast8Nz] / dnCols)*100.0;
    features[OFpercColsLeast16Nz] = (features[OFcolsLeast16Nz] / dnCols)*100.0;
    features[OFpercColsLeast32Nz] = (features[OFcolsLeast32Nz] / dnCols)*100.0;
    features[OFpercColsLeast64Nz] = (features[OFcolsLeast64Nz] / dnCols)*100.0;
    features[OFpercColsLeast128Nz] = (features[OFcolsLeast128Nz] / dnCols)*100.0;
    features[OFpercColsLeast256Nz] = (features[OFcolsLeast256Nz] / dnCols)*100.0;
    features[OFpercColsLeast512Nz] = (features[OFcolsLeast512Nz] / dnCols)*100.0;
    features[OFpercColsLeast1024Nz] = (features[OFcolsLeast1024Nz] / dnCols)*100.0;
    features[OFpercColsLeast2048Nz] = (features[OFcolsLeast2048Nz] / dnCols)*100.0;
    features[OFpercColsLeast4096Nz] = (features[OFcolsLeast4096Nz] / dnCols)*100.0;
}
