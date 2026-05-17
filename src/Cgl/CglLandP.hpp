// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     07/21/05
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------
#ifndef CglLandP_H
#define CglLandP_H

#include "CglLandPValidator.hpp"
#include "CglCutGenerator.hpp"
#include "CglParam.hpp"

#include <iostream>
class CoinWarmStartBasis;
/** Performs one round of Lift & Project using CglLandPSimplex
    to build cuts
*/

namespace LAP
{
enum LapMessagesTypes
{
    BEGIN_ROUND,
    END_ROUND,
    DURING_SEP,
    CUT_REJECTED,
    CUT_FAILED,
    CUT_GAP,
    LAP_CUT_FAILED_DO_MIG,
    LAP_MESSAGES_DUMMY_END
};
/** Output messages for Cgl */
class LapMessages : public CoinMessages
{
public:
    /** Constructor */
    LapMessages( );
    /** destructor.*/
    virtual ~LapMessages() {}
};
class CglLandPSimplex;
}

class CGLLIB_EXPORT CglLandP : public CglCutGenerator
{
    friend CGLLIB_EXPORT void CglLandPUnitTest(OsiSolverInterface *si, const std::string & mpsDir);

    friend class LAP::CglLandPSimplex;
    friend class CftCglp;

public:

    enum SelectionRules
    {
        mostNegativeRc /** select most negative reduced cost */,
        bestPivot /** select best possible pivot.*/,
        initialReducedCosts/** Select only those rows which had initialy a 0 reduced cost.*/
    };

    enum ExtraCutsMode
    {
        none/** Generate no extra cuts.*/,
        AtOptimalBasis /** Generate cuts from the optimal basis.*/,
        WhenEnteringBasis /** Generate cuts as soon as a structural enters the basis.*/,
        AllViolatedMigs/** Generate all violated Mixed integer Gomory cuts in the course of the optimization.*/
    };

    /** Space where cuts are optimized.*/
    enum SeparationSpaces
    {
        Fractional=0 /** True fractional space.*/,
        Fractional_rc/** Use fractional space only for computing reduced costs.*/,
        Full /** Work in full space.*/
    };

    /** Normalization */
    enum Normalization
    {
        Unweighted = 0,
        WeightRHS,
        WeightLHS,
        WeightBoth
    };

    enum LHSnorm
    {
        L1 = 0,
        L2,
        SupportSize,
        Infinity,
        Average,
        Uniform
    };
    /** RHS weight in normalization.*/
    enum RhsWeightType
    {
        Fixed = 0 /** 2*initial number of constraints. */,
        Dynamic /** 2 * current number of constraints. */
    };
    /** Class storing parameters.
        \remark I take all parameters from Ionut's code */
    class CGLLIB_EXPORT Parameters : public CglParam
    {
    public:
        /** Default constructor (with default values)*/
        Parameters();
        /** Copy constructor */
        Parameters(const Parameters &other);
        /** Assignment opertator */
        Parameters & operator=(const Parameters &other);
        /// @name integer parameters
        ///@{

        /** Max number of pivots before we generate the cut
          \default 20 */
        int pivotLimit;
        /** Max number of pivots at regular nodes. Put a value if you want it lower than the global pivot limit.
         \default 100.*/
        int pivotLimitInTree;
        /** Maximum number of cuts generated at a given round*/
        int maxCutPerRound;
        /** Maximum number of failed pivots before aborting */
        int failedPivotLimit;
        /** maximum number of consecutive degenerate pivots
          \default 0 */
        int degeneratePivotLimit;
        /** Maximum number of extra rows to generate per round.*/
        int extraCutsLimit;
	/// Maximum size of "indices"
	int maximumCandidates;
	/// Maximum size of cut
	int maximumCutLength;
        ///@}
        /// @name double parameters
        ///@{
        /** Tolerance for small pivots values (should be the same as the solver */
        double pivotTol;
        /** A variable have to be at least away from integrity to be generated */
        double away;
        /** Total time limit for cut generation.*/
        double timeLimit;
        /** Time limit for generating a single cut.*/
        double singleCutTimeLimit;
        /** Weight to put in RHS of normalization if static.*/
        double rhsWeight;
        ///@}

        /// @name Flags
        ///@{
        /** Do we use tableau row or the disjunction (I don't really get that there should be a way to always use the tableau)*/
        bool useTableauRow;
        /** Do we apply Egon Balas's Heuristic for modularized cuts */
        bool modularize;
        /** Do we strengthen the final cut (always do if modularize is 1) */
        bool strengthen;
        /** Wether to limit or not the number of mistaken RC (when perturbation is applied).*/
        bool countMistakenRc;
        /** Work in the reduced space (only non-structurals enter the basis) */
        SeparationSpaces sepSpace;
        /** Apply perturbation procedure. */
        bool perturb;
        /** How to weight normalization.*/
        Normalization normalization;
        /** How to weight RHS of normalization.*/
        RhsWeightType rhsWeightType;
        /** How to weight LHS of normalization.*/
        LHSnorm lhs_norm;
        /** Generate extra constraints from optimal lift-and-project basis.*/
        ExtraCutsMode generateExtraCuts;
        /** Which rule to apply for choosing entering and leaving variables.*/
        SelectionRules pivotSelection;
        ///@}
    };


    /** Constructor for the class*/
    CglLandP(const CglLandP::Parameters &params = CglLandP::Parameters(),
             const LAP::Validator &validator = LAP::Validator());
    /** Destructor */
    ~CglLandP();
    /** Copy constructor */
    CglLandP(const CglLandP &source);
    /** Assignment operator */
    CglLandP& operator=(const CglLandP &rhs);
    /** Clone function */
    CglCutGenerator * clone() const;

    /**@name Generate Cuts */
    //@{

    virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
                              const CglTreeInfo info = CglTreeInfo());

    //@}

    virtual bool needsOptimalBasis() const
    {
        return true;
    }

    LAP::Validator & validator()
    {
        return validator_;
    }
    /** set level of log for cut generation procedure :
    <ol start=0 >
    	<li> for none </li>
    	<li> for log at begin and end of procedure + at some time interval </li>
    	<li> for log at every cut generated </li>
    	</ol>
    	*/
    void setLogLevel(int level)
    {
        handler_->setLogLevel(level);
    }

    class NoBasisError : public CoinError
    {
    public:
        NoBasisError(): CoinError("No basis available","LandP","") {}
    };

    class SimplexInterfaceError : public CoinError
    {
    public:
        SimplexInterfaceError(): CoinError("Invalid conversion to simplex interface", "CglLandP","CglLandP") {}
    };
    Parameters & parameter()
    {
        return params_;
    }
private:


    void scanExtraCuts(OsiCuts& cs, const double * colsol) const;

    Parameters params_;

    /** Some informations that will be changed by the pivots and that we want to keep*/
    struct CachedData
    {
        CachedData(int nBasics = 0 , int nNonBasics = 0);
        CachedData(const CachedData & source);

        CachedData& operator=(const CachedData &source);
        /** Get the data from a problem */
        void getData(const OsiSolverInterface &si);

        void clean();

        ~CachedData();
        /** Indices of basic variables in starting basis (ordered if variable basics_[i] s basic in row i)*/
        int * basics_;
        /** Indices of non-basic variables */
        int *nonBasics_;
        /** number of basics variables */
        int nBasics_;
        /** number of non-basics */
        int nNonBasics_;
        /** Optimal basis */
        CoinWarmStartBasis * basis_;
        /** Stores the value of the solution to cut */
        double * colsol_;
        /** Stores the values of the slacks */
        double * slacks_;
        /** Stores wheter slacks are integer constrained */
        bool * integers_;
        /** Solver before pivots */
        OsiSolverInterface * solver_;
    };
    /** Retrieve sorted integer variables which are fractional in the solution.
        Return the number of variables.*/
    int getSortedFractionals(CoinPackedVector &xFrac,
                             const CachedData & data,
                             const CglLandP::Parameters& params) const;
    /** Retrieve sorted integer variables which are fractional in the solution.
        Return the number of variables.*/
    void getSortedFractionalIndices(std::vector<int>& indices,
                                    const CachedData &data,
                                    const CglLandP::Parameters & params) const;
    /** Cached informations about problem.*/
    CachedData cached_;
    /** message handler */
    CoinMessageHandler * handler_;
    /** messages */
    CoinMessages messages_;
    /** cut validator */
    LAP::Validator validator_;
    /** number of rows in the original problems. */
    int numrows_;
    /** number of columns in the original problems. */
    int numcols_;
    /** Original lower bounds for the problem (for lifting cuts).*/
    double * originalColLower_;
    /** Original upper bounds for the problem (for lifting cuts).*/
    double * originalColUpper_;
    /** Flag to say if cuts can be lifted.*/
    bool canLift_;
    /** Store some extra cut which could be cheaply generated but do not cut current incumbent.*/
    OsiCuts extraCuts_;
};
CGLLIB_EXPORT
void CglLandPUnitTest(OsiSolverInterface *si, const std::string & mpsDir);

#endif

