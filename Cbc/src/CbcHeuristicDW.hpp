// $Id: CbcHeuristicDW.hpp 1899 2013-04-09 18:12:08Z stefan $
// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).


#ifndef CbcHeuristicDW_H
#define CbcHeuristicDW_H

#include "CbcHeuristic.hpp"

/** 
    This is unlike the other heuristics in that it is very very compute intensive.
    It tries to find a DW structure and use that
 */

class CbcHeuristicDW : public CbcHeuristic {
public:

    // Default Constructor
    CbcHeuristicDW ();

    /* Constructor with model - assumed before cuts
    */
    CbcHeuristicDW (CbcModel & model, int keepContinuous=0);

    /* Constructor with model - assumed before cuts
    */
    CbcHeuristicDW (CbcModel & model,
		    int callBack(CbcHeuristicDW * currentHeuristic,
				 CbcModel * thisModel,
				 int whereFrom),
		    int keepContinuous=0);

    // Copy constructor
    CbcHeuristicDW ( const CbcHeuristicDW &);

    // Destructor
    ~CbcHeuristicDW ();

    /// Clone
    virtual CbcHeuristic * clone() const;


    /// Assignment operator
    CbcHeuristicDW & operator=(const CbcHeuristicDW& rhs);

    /// Create C++ lines to get to current state
    virtual void generateCpp( FILE * fp) ;

    /// Resets stuff if model changes
    virtual void resetModel(CbcModel * model);

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel * model);
    using CbcHeuristic::solution ;
    /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
        This does Relaxation Induced Neighborhood Search
    */
    virtual int solution(double & objectiveValue,
                         double * newSolution);
    /** Return number of blocks
      <=0 - no usable structure */
    inline int numberBlocks() const
    { return numberBlocks_;}
    /// Pass in a solution
    void passInSolution(const double * solution);
    /// Pass in continuous solution
    void passInContinuousSolution(const double * solution);
    /** DW Proposal actions
        fullDWEverySoOften -
        0 - off
        k - every k times solution gets better
     */
    void setProposalActions(int fullDWEverySoOften);
    /// Objective value when whichDw created
    double objectiveValueWhen(int whichDW) const;
    /// Number of columns in DW
    int numberColumnsDW(int whichDW) const;
    /// Solver
    inline OsiSolverInterface * solver() const
    { return solver_;}
    /// DW model (user must delete)
    OsiSolverInterface * DWModel(int whichDW) const;
    /// Best objective value
    inline double bestObjective() const
    { return bestObjective_;}
    /// Best solution found so far
    inline const double * bestSolution() const
    { return bestSolution_;}
    /// Continuous solution
    inline const double * continuousSolution() const
    { return continuousSolution_;}
    /// Reduced costs of fixed solution
    inline const double * fixedDj() const
    { return fixedDj_;}
    /// Objective at which DW updated
    inline const double * objectiveDW() const
    { return objectiveDW_;}
    /// Number of times we have added to DW model
    inline int numberDWTimes() const
    { return numberDWTimes_;}
    /// Number of columns in DW
    inline const int * numberColumnsDW() const
    { return numberColumnsDW_;}
    /// Set number of passes
    inline void setNumberPasses(int value)
    { numberPasses_ = value;}
    /// Set number of passes without better solution
    inline void setNumberBadPasses(int value)
    { numberBadPasses_ = value;}
    /// Set number free integers needed (Base value)
    inline void setNumberNeeded(int value)
    { nNeededBase_ = value;}
    /// Get number free integers needed (Base value)
    inline int getNumberNeeded() const
    {return nNeededBase_;}
    /// Set number free integers needed (Current value)
    inline void setCurrentNumberNeeded(int value)
    { nNeeded_ = value;}
    /// Get number free integers needed (Current value)
    inline int getCurrentNumberNeeded() const
    {return nNeeded_;}
    /// Set number nodes (could be done in callback) (Base value)
    inline void setNumberNodes(int value)
    { nNodesBase_ = value;}
    /// Get number nodes (could be done in callback) (Base value)
    inline int getNumberNodes() const
    {return nNodesBase_;}
    /// Set number nodes (could be done in callback) (Current value)
    inline void setCurrentNumberNodes(int value)
    { nNodes_ = value;}
    /// Get number nodes (could be done in callback) (Current value)
    inline int getCurrentNumberNodes() const
    {return nNodes_;}
    /// Set target objective
    inline void setTargetObjective(double value)
    { targetObjective_ = value;}
    /// Sets how often to do it
    inline void setHowOften(int value) {
        howOften_ = value;
    }
    /// Block for every row
    inline const int * whichRowBlock() const
    { return whichRowBlock_;}
    /// Block for every column
    inline const int * whichColumnBlock() const
    { return whichColumnBlock_;}
    /// Initial Lower bounds
    inline double * initialLower() const
    { return saveLower_;}
    /// Initial Upper bounds
    inline double * initialUpper() const
    { return saveUpper_;}
    /// Local integer arrays (each numberBlocks_ long)
    inline int * intArrays() const
    { return intArray_;}
    /// Local double arrays (each numberBlocks_ long)
    inline double * doubleArrays() const
    { return doubleArray_;}
    /// Phase of solution
    inline int phase() const
    { return phase_;}
    /// Pass number
    inline int pass() const
    { return pass_;}
    /// Which columns are in block
    inline const int * columnsInBlock() const
    { return columnsInBlock_;}
    /// Starts for columnsInBlock
    inline const int * startColumnBlock() const
    { return startColumnBlock_;}
    /// Number of integer variables in each block
    inline const int * intsInBlock() const
    { return intsInBlock_;}
    /// Objective value (could also check validity)
    double objectiveValue(const double * solution);
private:
    /// Guts of copy
    void gutsOfCopy(const CbcHeuristicDW & rhs);
    /// Guts of delete
    void gutsOfDelete();
    /// Set default values
    void setDefaults();
    /// Find structure
    void findStructure();
    /// Set up DW structure
    void setupDWStructures();
    /// Add DW proposals
    int addDW(const double * solution,int numberBlocksUsed, 
	      const int * whichBlocks);
protected:
    typedef int (*heuristicCallBack) (CbcHeuristicDW * ,CbcModel *, int) ;
    // Data
    /// Target objective
    double targetObjective_;
    /// Best objective value
    double bestObjective_;
    /// Objective value last time
    double lastObjective_;
    /** Call back
	whereFrom -
	0 - after blocks found but before data setup
	1 - after blocks sorted but before used
	2 - just before normal branch and bound
	3 - after DW has been updated
	4 - if better solution found
	5 - every time a block might be used
	next few for adjustment of nNeeded etc
	6 - complete search done - no solution
	7 - stopped on nodes - no improvement
	8 - improving (same as 4 but after nNeeded changed
	Pointers to local data given by following pointers
    */
    heuristicCallBack functionPointer_;
    /// Local integer arrays (each numberBlocks_ long)
    int * intArray_;
    /// Local double arrays (each numberBlocks_ long)
    double * doubleArray_;
    /// Base solver
    OsiSolverInterface * solver_;
    /// DW solver
    OsiSolverInterface * dwSolver_;
    /// Best solution found so far
    double * bestSolution_;
    /// Continuous solution
    double * continuousSolution_;
    /// Reduced costs of fixed solution
    double * fixedDj_;
    /// Original lower bounds
    double * saveLower_;
    /// Original Upper bounds
    double * saveUpper_;
    /// random numbers for master rows
    double * random_;
    /// Weights for each proposal
    double * weights_;
    /// Objective at which DW updated
    double * objectiveDW_;
    /// Number of columns in each DW
    int * numberColumnsDW_;
    /// Block for every row
    int * whichRowBlock_;
    /// Block for every column
    int * whichColumnBlock_;
    /// Block number for each proposal
    int * dwBlock_;
    /// Points back to master rows
    int * backwardRow_;
    /// Which rows are in blocke
    int * rowsInBlock_;
    /// Which columns are in block
    int * columnsInBlock_;
    /// Starts for rowsInBlock
    int * startRowBlock_;
    /// Starts for columnsInBlock
    int * startColumnBlock_;
    /// Number of integer variables in each block
    int * intsInBlock_;
    /// Bits set for 1 integers in each block 
    unsigned int * fingerPrint_;
    /// Affinity each block has for other (will be triangular?)
    unsigned short * affinity_;
    /** DW Proposal actions
        fullDWEverySoOften -
        0 - off
        k - every k times solution gets better
     */
    int fullDWEverySoOften_;
    /// Number of passes
    int numberPasses_;
    /// How often to do (code can change)
    int howOften_;
    /// Current maximum number of DW proposals
    int maximumDW_;
    /// Number of DW proposals
    int numberDW_;
    /// Number of times we have added to DW model
    int numberDWTimes_;
    /// Number of unsigned ints needed for each block of fingerPrint
    int sizeFingerPrint_;
    /// Number of columns in master
    int numberMasterColumns_;
    /// Number of rows in master
    int numberMasterRows_;
    /// Number of blocks
    int numberBlocks_;
    /// Action on decomposition - 1 keep continuous, 0 don't
    int keepContinuous_;
    /// Phase of solution
    int phase_;
    /// Pass number
    int pass_;
    /// Base number of integers needed
    int nNeededBase_;
    /// Base number of nodes needed
    int nNodesBase_;
    /// Base number of integers needed
    int nNeeded_;
    /// Base number of nodes needed
    int nNodes_;
    /// Number of passes without better solution
    int numberBadPasses_;
    // 0 - fine, 1 can't be better, 2 max node
    int solveState_;
};

#endif
