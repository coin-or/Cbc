// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
/* $Id: ClpAmplStuff.cpp 1200 2009-07-25 08:44:13Z forrest $ */

#include "ClpConfig.h"
#include "CbcConfig.h"
#ifdef COIN_HAS_ASL
#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpSimplex.hpp"
#include "ClpAmplObjective.hpp"
#include "ClpConstraintAmpl.hpp"
#include "ClpMessage.hpp"
#include "CoinUtilsConfig.h"
#include "CoinHelperFunctions.hpp"
#include "CoinWarmStartBasis.hpp"
#include "OsiSolverInterface.hpp"
#include "CbcSolver.hpp"
#include "Cbc_ampl.h"
#include "CoinTime.hpp"
#include "CglStored.hpp"
#include "CoinModel.hpp"
#include "CbcLinked.hpp"
class CbcAmpl  : public CbcUser {

public:
    ///@name usage methods
    //@{
    /// Solve (whatever that means)
    virtual void solve(CbcSolver * model, const char * options);
    /// Returns true if function knows about option
    virtual bool canDo(const char * options) ;
    /** Import - gets full command arguments
        Returns -1 - no action
                 0 - data read in without error
             1 - errors
    */
    virtual int importData(CbcSolver * model, int & argc, char ** & argv);
    /// Export 1 OsiClpSolver, 2 CbcModel - add 10 if infeasible from odd situation
    virtual void exportSolution(CbcSolver * model, int mode, const char * message = NULL) ;
    /// Export Data (i.e. at very end)
    virtual void exportData(CbcSolver * model);
    /// Get useful stuff
    virtual void fillInformation(CbcSolver * model,
                                 CbcSolverUsefulData & info);
    //@}
    ///@name Constructors and destructors etc
    //@{
    /// Default Constructor
    CbcAmpl();

    /** Copy constructor .
     */
    CbcAmpl(const CbcAmpl & rhs);

    /// Assignment operator
    CbcAmpl & operator=(const CbcAmpl& rhs);

    /// Clone
    virtual CbcUser * clone() const;

    /// Destructor
    virtual ~CbcAmpl ();
    //@}
private:
    ///@name Private member data
    //@{
    /// AMPL info
    ampl_info info_;
    //@}
};
// Mechanicsburg stuff
CbcAmpl::CbcAmpl()
        : CbcUser()
{
    userName_ = "mech";
    memset(&info_, 0, sizeof(info_));
}
CbcAmpl::~CbcAmpl()
{
}
// Copy constructor
CbcAmpl::CbcAmpl ( const CbcAmpl & rhs)
        : CbcUser(rhs)
{
    info_ = rhs.info_;
}
// Assignment operator
CbcAmpl &
CbcAmpl::operator=(const CbcAmpl & rhs)
{
    if (this != &rhs) {
        CbcUser::operator=(rhs);
        info_ = rhs.info_;
    }
    return *this;
}
// Clone
CbcUser *
CbcAmpl::clone() const
{
    return new CbcAmpl(*this);
}
// Solve (whatever that means)
void
CbcAmpl::solve(CbcSolver * controlModel, const char * options)
{
    assert (controlModel->model());
    //OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model->solver());
    //ClpSimplex * lpSolver = clpSolver->getModelPtr();
    if (!strcmp(options, "cbc_load")) {
    } else if (!strcmp(options, "cbc_quit")) {
    } else {
        printf("unknown option for CbcAmpl is %s\n", options);
        abort();
    }
}
// Returns true if function knows about option
bool
CbcAmpl::canDo(const char * options)
{
    return (!strcmp(options, "cbc_load") || !strcmp(options, "cbc_quit"));
}
/* Import - gets full command arguments
   Returns -1 - no action
            0 - data read in without error
	    1 - errors
*/
int
CbcAmpl::importData(CbcSolver * control, int &argc, char ** & argv)
{
    CbcModel * babModel = control->model();
    assert (babModel);
    CoinMessageHandler * generalMessageHandler = babModel->messageHandler();
    OsiClpSolverInterface * solver = dynamic_cast< OsiClpSolverInterface*> (control->model()->solver());
    assert (solver);
    CoinMessages generalMessages = solver->getModelPtr()->messages();
    char generalPrint[10000];
    OsiSolverLink * si = NULL;
    ClpSimplex * lpSolver = solver->getModelPtr();
    if (argc > 2 && !strcmp(argv[2], "-AMPL")) {
        // see if log in list
        bool printing = false;
        for (int i = 1; i < argc; i++) {
            if (!strncmp(argv[i], "log", 3)) {
                const char * equals = strchr(argv[i], '=');
                if (equals && atoi(equals + 1) > 0) {
                    printing = true;
                    info_.logLevel = atoi(equals + 1);
                    control->setIntValue(LOGLEVEL, info_.logLevel);
                    // mark so won't be overWritten
                    info_.numberRows = -1234567;
                    break;
                }
            }
        }
        union {
            void * voidModel;
            CoinModel * model;
        } coinModelStart;
        coinModelStart.model = NULL;
        int returnCode = readAmpl(&info_, argc, argv, & coinModelStart.voidModel);
        CoinModel * coinModel = coinModelStart.model;
        if (returnCode)
            return returnCode;
        control->setReadMode(3); // so will start with parameters
        // see if log in list (including environment)
        for (int i = 1; i < info_.numberArguments; i++) {
            if (!strcmp(info_.arguments[i], "log")) {
                if (i < info_.numberArguments - 1 && atoi(info_.arguments[i+1]) > 0)
                    printing = true;
                break;
            }
        }
        control->setPrinting(printing);
        if (printing)
            printf("%d rows, %d columns and %d elements\n",
                   info_.numberRows, info_.numberColumns, info_.numberElements);
        if (!coinModel) {
            solver->loadProblem(info_.numberColumns, info_.numberRows, info_.starts,
                                info_.rows, info_.elements,
                                info_.columnLower, info_.columnUpper, info_.objective,
                                info_.rowLower, info_.rowUpper);
            if (info_.numberSos) {
                // SOS
                solver->setSOSData(info_.numberSos, info_.sosType, info_.sosStart,
                                   info_.sosIndices, info_.sosReference);
            }
        } else {
            // save
            control->setOriginalCoinModel(coinModel);
            // load from coin model
            OsiSolverLink solver1;
            OsiSolverInterface * solver2 = solver1.clone();
            babModel->assignSolver(solver2, false);
            si = dynamic_cast<OsiSolverLink *>(babModel->solver()) ;
            assert (si != NULL);
            si->setDefaultMeshSize(0.001);
            // need some relative granularity
            si->setDefaultBound(100.0);
            double dextra3 = control->doubleValue(DEXTRA3);
            if (dextra3)
                si->setDefaultMeshSize(dextra3);
            si->setDefaultBound(100000.0);
            si->setIntegerPriority(1000);
            si->setBiLinearPriority(10000);
            CoinModel * model2 = (CoinModel *) coinModel;
            int logLevel = control->intValue(LOGLEVEL);
            si->load(*model2, true, logLevel);
            // redo
            solver = dynamic_cast< OsiClpSolverInterface*> (control->model()->solver());
            lpSolver = solver->getModelPtr();
            solver->messageHandler()->setLogLevel(0) ;
            control->setIntValue(TESTOSI, 0);
            if (info_.cut) {
                printf("Sorry - can't do cuts with LOS as ruins delicate row order\n");
                abort();
            }
        }
        if (info_.cut) {
            int numberRows = info_.numberRows;
            int * whichRow = new int [numberRows];
            // Row copy
            const CoinPackedMatrix * matrixByRow = solver->getMatrixByRow();
            const double * elementByRow = matrixByRow->getElements();
            const int * column = matrixByRow->getIndices();
            const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
            const int * rowLength = matrixByRow->getVectorLengths();

            const double * rowLower = solver->getRowLower();
            const double * rowUpper = solver->getRowUpper();
            int nDelete = 0;
            CglStored storedAmpl;
            for (int iRow = 0; iRow < numberRows; iRow++) {
                if (info_.cut[iRow]) {
                    whichRow[nDelete++] = iRow;
                    int start = rowStart[iRow];
                    storedAmpl.addCut(rowLower[iRow], rowUpper[iRow],
                                      rowLength[iRow], column + start, elementByRow + start);
                }
            }
            control->addCutGenerator(&storedAmpl);
            solver->deleteRows(nDelete, whichRow);
            // and special matrix
            si->cleanMatrix()->deleteRows(nDelete, whichRow);
            delete [] whichRow;
        }
        // If we had a solution use it
        if (info_.primalSolution) {
            solver->setColSolution(info_.primalSolution);
        }
        // status
        if (info_.rowStatus) {
            unsigned char * statusArray = lpSolver->statusArray();
            memset(statusArray, 0, lpSolver->numberColumns() + lpSolver->numberRows());
            int i;
            for (i = 0; i < info_.numberColumns; i++)
                statusArray[i] = (char)info_.columnStatus[i];
            statusArray += info_.numberColumns;
            for (i = 0; i < info_.numberRows; i++)
                statusArray[i] = (char)info_.rowStatus[i];
            CoinWarmStartBasis * basis = lpSolver->getBasis();
            solver->setWarmStart(basis);
            delete basis;
        }
        freeArrays1(&info_);
        // modify objective if necessary
        solver->setObjSense(info_.direction);
        solver->setDblParam(OsiObjOffset, info_.offset);
        if (info_.offset) {
            sprintf(generalPrint, "Ampl objective offset is %g",
                    info_.offset);
            generalMessageHandler->message(CLP_GENERAL, generalMessages)
            << generalPrint
            << CoinMessageEol;
        }
        // Set integer variables (unless nonlinear when set)
        if (!info_.nonLinear) {
            for (int i = info_.numberColumns - info_.numberIntegers;
                    i < info_.numberColumns; i++)
                solver->setInteger(i);
        }
        // change argc etc
        argc = info_.numberArguments;
        argv = info_.arguments;
        return 0;
    } else {
        return -1;
    }
    abort();
    return -1;
}
// Export 1 OsiClpSolver, 2 CbcModel - add 10 if infeasible from odd situation
void
CbcAmpl::exportSolution(CbcSolver * model, int mode, const char * message)
{
    OsiClpSolverInterface * solver = model->originalSolver();
    if (!solver) {
        solver = dynamic_cast< OsiClpSolverInterface*> (model->model()->solver());
        assert (solver);
    }
    ClpSimplex * lpSolver = solver->getModelPtr();
    int numberColumns = lpSolver->numberColumns();
    int numberRows = lpSolver->numberRows();
    double totalTime = CoinCpuTime() - model->startTime();
    if (mode == 1) {
        double value = lpSolver->getObjValue() * lpSolver->getObjSense();
        char buf[300];
        int pos = 0;
        int iStat = lpSolver->status();
        if (iStat == 0) {
            pos += sprintf(buf + pos, "optimal," );
        } else if (iStat == 1) {
            // infeasible
            pos += sprintf(buf + pos, "infeasible,");
        } else if (iStat == 2) {
            // unbounded
            pos += sprintf(buf + pos, "unbounded,");
        } else if (iStat == 3) {
            pos += sprintf(buf + pos, "stopped on iterations or time,");
        } else if (iStat == 4) {
            iStat = 7;
            pos += sprintf(buf + pos, "stopped on difficulties,");
        } else if (iStat == 5) {
            iStat = 3;
            pos += sprintf(buf + pos, "stopped on ctrl-c,");
        } else {
            pos += sprintf(buf + pos, "status unknown,");
            iStat = 6;
        }
        info_.problemStatus = iStat;
        info_.objValue = value;
        pos += sprintf(buf + pos, " objective %.*g", ampl_obj_prec(),
                       value);
        sprintf(buf + pos, "\n%d iterations",
                lpSolver->getIterationCount());
        free(info_.primalSolution);
        info_.primalSolution = (double *) malloc(numberColumns * sizeof(double));
        CoinCopyN(lpSolver->primalColumnSolution(), numberColumns, info_.primalSolution);
        free(info_.dualSolution);
        info_.dualSolution = (double *) malloc(numberRows * sizeof(double));
        CoinCopyN(lpSolver->dualRowSolution(), numberRows, info_.dualSolution);
        CoinWarmStartBasis * basis = lpSolver->getBasis();
        free(info_.rowStatus);
        info_.rowStatus = (int *) malloc(numberRows * sizeof(int));
        free(info_.columnStatus);
        info_.columnStatus = (int *) malloc(numberColumns * sizeof(int));
        // Put basis in
        int i;
        // free,basic,ub,lb are 0,1,2,3
        for (i = 0; i < numberRows; i++) {
            CoinWarmStartBasis::Status status = basis->getArtifStatus(i);
            info_.rowStatus[i] = status;
        }
        for (i = 0; i < numberColumns; i++) {
            CoinWarmStartBasis::Status status = basis->getStructStatus(i);
            info_.columnStatus[i] = status;
        }
        // put buffer into info_
        strcpy(info_.buffer, buf);
        delete basis;
    } else if (mode == 2) {
        CbcModel * babModel = model->model();
        int iStat = babModel->status();
        int iStat2 = babModel->secondaryStatus();
        double value = babModel->getObjValue() * lpSolver->getObjSense();
        char buf[300];
        int pos = 0;
        if (iStat == 0) {
            if (babModel->getObjValue() < 1.0e40) {
                pos += sprintf(buf + pos, "optimal," );
            } else {
                // infeasible
                iStat = 1;
                pos += sprintf(buf + pos, "infeasible,");
            }
        } else if (iStat == 1) {
            if (iStat2 != 6)
                iStat = 3;
            else
                iStat = 4;
            std::string minor[] = {"", "", "gap", "nodes", "time", "", "solutions", "user ctrl-c"};
            pos += sprintf(buf + pos, "stopped on %s,", minor[iStat2].c_str());
        } else if (iStat == 2) {
            iStat = 7;
            pos += sprintf(buf + pos, "stopped on difficulties,");
        } else if (iStat == 5) {
            iStat = 3;
            pos += sprintf(buf + pos, "stopped on ctrl-c,");
        } else {
            pos += sprintf(buf + pos, "status unknown,");
            iStat = 6;
        }
        info_.problemStatus = iStat;
        info_.objValue = value;
        if (babModel->getObjValue() < 1.0e40) {
            int precision = ampl_obj_prec();
            if (precision > 0)
                pos += sprintf(buf + pos, " objective %.*g", precision,
                               value);
            else
                pos += sprintf(buf + pos, " objective %g", value);
        }
        sprintf(buf + pos, "\n%d nodes, %d iterations, %g seconds",
                babModel->getNodeCount(),
                babModel->getIterationCount(),
                totalTime);
        if (babModel->bestSolution()) {
            free(info_.primalSolution);
            info_.primalSolution = (double *) malloc(numberColumns * sizeof(double));
            CoinCopyN(lpSolver->primalColumnSolution(), numberColumns, info_.primalSolution);
            free(info_.dualSolution);
            info_.dualSolution = (double *) malloc(numberRows * sizeof(double));
            CoinCopyN(lpSolver->dualRowSolution(), numberRows, info_.dualSolution);
        } else {
            info_.primalSolution = NULL;
            info_.dualSolution = NULL;
        }
        // put buffer into info
        strcpy(info_.buffer, buf);
    } else if (mode == 11 || mode == 12) {
        // infeasible
        info_.problemStatus = 1;
        info_.objValue = 1.0e100;
        sprintf(info_.buffer, "%s", message);
        info_.primalSolution = NULL;
        info_.dualSolution = NULL;
    }
}
// Export Data (i.e. at very end)
void
CbcAmpl::exportData(CbcSolver * model)
{
    writeAmpl(&info_);
    freeArrays2(&info_);
    freeArgs(&info_);
}
// Get useful stuff
void
CbcAmpl::fillInformation(CbcSolver * model,
                         CbcSolverUsefulData & info)
{
    memset(&info, 0, sizeof(info));
    info.priorities_ = info_.priorities;
    info.sosPriority_ = info_.sosPriority;
    info.branchDirection_ = info_.branchDirection;
    info.primalSolution_ = info_.primalSolution;
    info.pseudoDown_ = info_.pseudoDown;
    info.pseudoUp_ = info_.pseudoUp;
}
void addAmplToCbc(CbcSolver * control)
{
    CbcAmpl ampl;
    control->addUserFunction(&ampl);
}
extern "C" {
    //# include "getstub.h"
# include "asl_pfgh.h"
}
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpAmplObjective::ClpAmplObjective ()
        : ClpObjective()
{
    type_ = 12;
    objective_ = NULL;
    amplObjective_ = NULL;
    gradient_ = NULL;
    offset_ = 0.0;
}
// stolen from IPopt with changes
typedef struct {
    double obj_sign_;
    ASL_pfgh * asl_;
    double * non_const_x_;
    int * column_; // for jacobian
    int * rowStart_;
    double * gradient_;
    double * constraintValues_;
    int nz_h_full_; // number of nonzeros in hessian
    int nerror_;
    bool objval_called_with_current_x_;
    bool conval_called_with_current_x_;
    bool jacval_called_with_current_x_;
} CbcAmplInfo;
#if 0
static bool get_nlp_info(void * amplInfo, int & n, int & m, int & nnz_jac_g,
                         int & nnz_h_lag)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;

    n = n_var; // # of variables
    m = n_con; // # of constraints
    nnz_jac_g = nzc; // # of non-zeros in the jacobian
    nnz_h_lag = info->nz_h_full_; // # of non-zeros in the hessian

    return true;
}

static bool get_bounds_info(void * amplInfo, int  n, double * x_l,
                            double * x_u, int  m, double * g_l, double * g_u)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    assert(n == n_var);
    assert(m == n_con);
    int i;
    for (i = 0; i < n; i++) {
        x_l[i] = LUv[2*i];
        x_u[i] = LUv[2*i+1];
    }

    for (i = 0; i < m; i++) {
        g_l[i] = LUrhs[2*i];
        g_u[i] = LUrhs[2*i+1];
    }
    return true;
}

#endif
bool get_constraints_linearity(void * amplInfo, int  n,
                               int * const_types)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    //check that n is good
    assert(n == n_con);
    // check that there are no network constraints
    assert(nlnc == 0 && lnc == 0);
    //the first nlc constraints are non linear the rest is linear
    int i;
    for (i = 0; i < nlc; i++) {
        const_types[i] = 1;
    }
    // the rest is linear
    for (i = nlc; i < n_con; i++)
        const_types[i] = 0;
    return true;
}
#if 0
bool get_starting_point(int  n, bool init_x, double * x, bool init_z,
                        double * z_L, double * z_U, int  m, bool init_lambda, double * lambda)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    assert(n == n_var);
    assert(m == n_con);
    int i;

    if (init_x) {
        for (i = 0; i < n; i++) {
            if (havex0[i]) {
                x[i] = X0[i];
            } else {
                x[i] = 0.0;
            }
        }
    }

    if (init_z) {
        for (i = 0; i < n; i++) {
            z_L[i] = z_U[i] = 1.0;
        }
    }

    if (init_lambda) {
        for (i = 0; i < m; i++) {
            if (havepi0[i]) {
                lambda[i] = pi0[i];
            } else {
                lambda[i] = 0.0;
            }
        }
    }

    return true;
}
#endif
static bool internal_objval(CbcAmplInfo * info , double & obj_val)
{
    ASL_pfgh* asl = info->asl_;
    info->objval_called_with_current_x_ = false; // in case the call below fails

    if (n_obj == 0) {
        obj_val = 0;
        info->objval_called_with_current_x_ = true;
        return true;
    }  else {
        double  retval = objval(0, info->non_const_x_, (fint*)info->nerror_);
        if (!info->nerror_) {
            obj_val = info->obj_sign_ * retval;
            info->objval_called_with_current_x_ = true;
            return true;
        } else {
            abort();
        }
    }

    return false;
}
static bool internal_conval(CbcAmplInfo * info , double * g)
{
    ASL_pfgh* asl = info->asl_;
    info->conval_called_with_current_x_ = false; // in case the call below fails
    assert (g);

    conval(info->non_const_x_, g, (fint*)info->nerror_);

    if (!info->nerror_) {
        info->conval_called_with_current_x_ = true;
        return true;
    } else {
        abort();
    }
    return false;
}

static bool apply_new_x(CbcAmplInfo * info  , bool new_x, int  n, const double * x)
{
    ASL_pfgh* asl = info->asl_;

    if (new_x) {
        // update the flags so these methods are called
        // before evaluating the hessian
        info->conval_called_with_current_x_ = false;
        info->objval_called_with_current_x_ = false;
        info->jacval_called_with_current_x_ = false;

        //copy the data to the non_const_x_
        if (!info->non_const_x_) {
            info->non_const_x_ = new double [n];
        }

        for (int  i = 0; i < n; i++) {
            info->non_const_x_[i] = x[i];
        }

        // tell ampl that we have a new x
        xknowne(info->non_const_x_, (fint*)info->nerror_);
        return info->nerror_ ? false : true;
    }

    return true;
}
static bool eval_f(void * amplInfo, int  n, const double * x, bool new_x, double & obj_value)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    if (!apply_new_x(info, new_x, n, x)) {
        return false;
    }

    return internal_objval(info, obj_value);
}

static bool eval_grad_f(void * amplInfo, int  n, const double * x, bool new_x, double * grad_f)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    if (!apply_new_x(info, new_x, n, x)) {
        return false;
    }
    int i;

    if (n_obj == 0) {
        for (i = 0; i < n; i++) {
            grad_f[i] = 0.;
        }
    } else {
        objgrd(0, info->non_const_x_, grad_f, (fint*)info->nerror_);
        if (info->nerror_) {
            return false;
        }

        if (info->obj_sign_ == -1) {
            for (i = 0; i < n; i++) {
                grad_f[i] = -grad_f[i];
            }
        }
    }
    return true;
}
static bool eval_g(void * amplInfo, int  n, const double * x, bool new_x, double * g)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
#ifndef NDEBUG
    ASL_pfgh* asl = info->asl_;
#endif
    // warning: n_var is a macro that assumes we have a variable called asl
    assert(n == n_var);

    if (!apply_new_x(info, new_x, n, x)) {
        return false;
    }

    return internal_conval(info, g);
}

static bool eval_jac_g(void * amplInfo, int  n, const double * x, bool new_x,
                       double * values)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    assert(n == n_var);

    assert (values);
    if (!apply_new_x(info, new_x, n, x)) {
        return false;
    }

    jacval(info->non_const_x_, values, (fint*)info->nerror_);
    if (!info->nerror_) {
        return true;
    } else {
        abort();
    }
    return false;
}
#if 0

static bool eval_h(void * amplInfo, int  n, const double * x, bool new_x,
                   double  obj_factor, int  m, const double * lambda,
                   bool new_lambda, int  nele_hess, int * iRow,
                   int * jCol, double * values)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    assert(n == n_var);
    assert(m == n_con);
    int i;

    if (iRow && jCol && !values) {
        // setup the structure
        int k = 0;
        for (int i = 0; i < n; i++) {
            for (int j = sputinfo->hcolstarts[i]; j < sputinfo->hcolstarts[i+1]; j++) {
                //iRow[k] = i + 1;
                iRow[k] = i;
                jCol[k] = sputinfo->hrownos[j] + 1;
                k++;
            }
        }
        assert(k == nele_hess);
        return true;
    } else if (!iRow & !jCol && values) {
        if (!apply_new_x(info, new_x, n, x)) {
            return false;
        }
        if (!info->objval_called_with_current_x_) {
            double  dummy;
            internal_objval(info, dummy);
            internal_conval(info, m, NULL);
        }
        if (!info->conval_called_with_current_x_) {
            internal_conval(info, m, NULL);
        }
        // copy lambda to non_const_lambda - note, we do not store a copy like
        // we do with x since lambda is only used here and not in other calls
        double * non_const_lambda = new double [m];
        for (i = 0; i < m; i++) {
            non_const_lambda[i] = lambda[i];
        }

        real ow = info->obj_sign_ * obj_factor;
        sphes(values, -1, &ow, non_const_lambda);

        delete [] non_const_lambda;
        return true;
    } else {
        assert(false && "Invalid combination of iRow, jCol, and values pointers");
    }

    return false;
}
#endif
//-------------------------------------------------------------------
// Useful Constructor
//-------------------------------------------------------------------
ClpAmplObjective::ClpAmplObjective (void * amplInfo)
        : ClpObjective()
{
    type_ = 12;
    activated_ = 1;
    gradient_ = NULL;
    objective_ = NULL;
    offset_ = 0.0;
    amplObjective_ = amplInfo;
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpAmplObjective::ClpAmplObjective (const ClpAmplObjective & rhs)
        : ClpObjective(rhs)
{
    amplObjective_ = rhs.amplObjective_;
    offset_ = rhs.offset_;
    type_ = rhs.type_;
    if (!amplObjective_) {
        objective_ = NULL;
        gradient_ = NULL;
    } else {
        CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
        ASL_pfgh* asl = info->asl_;

        int numberColumns = n_var;;
        if (rhs.objective_) {
            objective_ = new double [numberColumns];
            memcpy(objective_, rhs.objective_, numberColumns*sizeof(double));
        } else {
            objective_ = NULL;
        }
        if (rhs.gradient_) {
            gradient_ = new double [numberColumns];
            memcpy(gradient_, rhs.gradient_, numberColumns*sizeof(double));
        } else {
            gradient_ = NULL;
        }
    }
}


//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpAmplObjective::~ClpAmplObjective ()
{
    delete [] objective_;
    delete [] gradient_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpAmplObjective &
ClpAmplObjective::operator=(const ClpAmplObjective & rhs)
{
    if (this != &rhs) {
        delete [] objective_;
        delete [] gradient_;
        amplObjective_ = rhs.amplObjective_;
        offset_ = rhs.offset_;
        type_ = rhs.type_;
        if (!amplObjective_) {
            objective_ = NULL;
            gradient_ = NULL;
        } else {
            CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
            ASL_pfgh* asl = info->asl_;

            int numberColumns = n_var;;
            if (rhs.objective_) {
                objective_ = new double [numberColumns];
                memcpy(objective_, rhs.objective_, numberColumns*sizeof(double));
            } else {
                objective_ = NULL;
            }
            if (rhs.gradient_) {
                gradient_ = new double [numberColumns];
                memcpy(gradient_, rhs.gradient_, numberColumns*sizeof(double));
            } else {
                gradient_ = NULL;
            }
        }
    }
    return *this;
}

// Returns gradient
double *
ClpAmplObjective::gradient(const ClpSimplex * model,
                           const double * solution, double & offset, bool refresh,
                           int includeLinear)
{
    if (model)
        assert (model->optimizationDirection() == 1.0);
    bool scaling = false;
    if (model && (model->rowScale() ||
                  model->objectiveScale() != 1.0 || model->optimizationDirection() != 1.0))
        scaling = true;
    const double * cost = NULL;
    if (model)
        cost = model->costRegion();
    if (!cost) {
        // not in solve
        cost = objective_;
        scaling = false;
    }
    assert (!scaling);
    if (!amplObjective_ || !solution || !activated_) {
        offset = offset_;
        return objective_;
    } else {
        if (refresh || !gradient_) {
            CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
            ASL_pfgh* asl = info->asl_;
            int numberColumns = n_var;;

            if (!gradient_)
                gradient_ = new double[numberColumns];
            assert (solution);
            eval_grad_f(amplObjective_, numberColumns, solution, true, gradient_);
            // Is this best way?
            double objValue = 0.0;
            eval_f(amplObjective_, numberColumns, solution, false, objValue);
            double objValue2 = 0.0;
            for (int i = 0; i < numberColumns; i++)
                objValue2 += gradient_[i] * solution[i];
            offset_ = objValue2 - objValue; // or other way???
            if (model && model->optimizationDirection() != 1.0) {
                offset *= model->optimizationDirection();
                for (int i = 0; i < numberColumns; i++)
                    gradient_[i] *= -1.0;
            }
        }
        offset = offset_;
        return gradient_;
    }
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpObjective * ClpAmplObjective::clone() const
{
    return new ClpAmplObjective(*this);
}
// Resize objective
void
ClpAmplObjective::resize(int newNumberColumns)
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    ASL_pfgh* asl = info->asl_;
    int numberColumns = n_var;;
    if (numberColumns != newNumberColumns) {
        abort();
    }

}
// Delete columns in  objective
void
ClpAmplObjective::deleteSome(int numberToDelete, const int * which)
{
    if (numberToDelete)
        abort();
}
/* Returns reduced gradient.Returns an offset (to be added to current one).
 */
double
ClpAmplObjective::reducedGradient(ClpSimplex * model, double * region,
                                  bool useFeasibleCosts)
{
    int numberRows = model->numberRows();
    int numberColumns = model->numberColumns();

    //work space
    CoinIndexedVector  * workSpace = model->rowArray(0);

    CoinIndexedVector arrayVector;
    arrayVector.reserve(numberRows + 1);

    int iRow;
#ifdef CLP_DEBUG
    workSpace->checkClear();
#endif
    double * array = arrayVector.denseVector();
    int * index = arrayVector.getIndices();
    int number = 0;
    const double * costNow = gradient(model, model->solutionRegion(), offset_,
                                      true, useFeasibleCosts ? 2 : 1);
    double * cost = model->costRegion();
    const int * pivotVariable = model->pivotVariable();
    for (iRow = 0; iRow < numberRows; iRow++) {
        int iPivot = pivotVariable[iRow];
        double value;
        if (iPivot < numberColumns)
            value = costNow[iPivot];
        else if (!useFeasibleCosts)
            value = cost[iPivot];
        else
            value = 0.0;
        if (value) {
            array[iRow] = value;
            index[number++] = iRow;
        }
    }
    arrayVector.setNumElements(number);

    // Btran basic costs
    model->factorization()->updateColumnTranspose(workSpace, &arrayVector);
    double * work = workSpace->denseVector();
    ClpFillN(work, numberRows, 0.0);
    // now look at dual solution
    double * rowReducedCost = region + numberColumns;
    double * dual = rowReducedCost;
    const double * rowCost = cost + numberColumns;
    for (iRow = 0; iRow < numberRows; iRow++) {
        dual[iRow] = array[iRow];
    }
    double * dj = region;
    ClpDisjointCopyN(costNow, numberColumns, dj);

    model->transposeTimes(-1.0, dual, dj);
    for (iRow = 0; iRow < numberRows; iRow++) {
        // slack
        double value = dual[iRow];
        value += rowCost[iRow];
        rowReducedCost[iRow] = value;
    }
    return offset_;
}
/* Returns step length which gives minimum of objective for
   solution + theta * change vector up to maximum theta.

   arrays are numberColumns+numberRows
*/
double
ClpAmplObjective::stepLength(ClpSimplex * model,
                             const double * solution,
                             const double * change,
                             double maximumTheta,
                             double & currentObj,
                             double & predictedObj,
                             double & thetaObj)
{
    // Assume convex
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    ASL_pfgh* asl = info->asl_;

    int numberColumns = n_var;;
    double * tempSolution = new double [numberColumns];
    double * tempGradient = new double [numberColumns];
    // current
    eval_f(amplObjective_, numberColumns, solution, true, currentObj);
    double objA = currentObj;
    double thetaA = 0.0;
    // at maximum
    int i;
    for (i = 0; i < numberColumns; i++)
        tempSolution[i] = solution[i] + maximumTheta * change[i];
    eval_f(amplObjective_, numberColumns, tempSolution, true, thetaObj);
    double objC = thetaObj;
    double thetaC = maximumTheta;
    double objB = 0.5 * (objA + objC);
    double thetaB = 0.5 * maximumTheta;
    double gradientNorm = 1.0e6;
    while (gradientNorm > 1.0e-6 && thetaC - thetaA > 1.0e-8) {
        for (i = 0; i < numberColumns; i++)
            tempSolution[i] = solution[i] + thetaB * change[i];
        eval_grad_f(amplObjective_, numberColumns, tempSolution, true, tempGradient);
        eval_f(amplObjective_, numberColumns, tempSolution, false, objB);
        double changeObj = 0.0;
        gradientNorm = 0.0;
        for (i = 0; i < numberColumns; i++) {
            changeObj += tempGradient[i] * change[i];
            gradientNorm += tempGradient[i] * tempGradient[i];
        }
        gradientNorm = fabs(changeObj) / sqrt(gradientNorm);
        // Should try and get quadratic convergence by interpolation
        if (changeObj < 0.0) {
            // increasing is good
            thetaA = thetaB;
        } else {
            // decreasing is good
            thetaC = thetaB;
        }
        thetaB = 0.5 * (thetaA + thetaC);
    }
    delete [] tempSolution;
    delete [] tempGradient;
    predictedObj = objB;
    return thetaB;
}
// Return objective value (without any ClpModel offset) (model may be NULL)
double
ClpAmplObjective::objectiveValue(const ClpSimplex * model, const double * solution) const
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    ASL_pfgh* asl = info->asl_;

    int numberColumns = n_var;;
    // current
    double currentObj = 0.0;
    eval_f(amplObjective_, numberColumns, solution, true, currentObj);
    return currentObj;
}
// Scale objective
void
ClpAmplObjective::reallyScale(const double * columnScale)
{
    abort();
}
/* Given a zeroed array sets nonlinear columns to 1.
   Returns number of nonlinear columns
*/
int
ClpAmplObjective::markNonlinear(char * which)
{
    int iColumn;
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    ASL_pfgh* asl = info->asl_;
    int nonLinear = CoinMax(nlvc, nlvo);
    for (iColumn = 0; iColumn < nonLinear; iColumn++) {
        which[iColumn] = 1;
    }
    int numberNonLinearColumns = 0;
    int numberColumns = n_var;;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (which[iColumn])
            numberNonLinearColumns++;
    }
    return numberNonLinearColumns;
}
// Say we have new primal solution - so may need to recompute
void
ClpAmplObjective::newXValues()
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplObjective_;
    info->conval_called_with_current_x_ = false;
    info->objval_called_with_current_x_ = false;
    info->jacval_called_with_current_x_ = false;
}
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################
//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpConstraintAmpl::ClpConstraintAmpl ()
        : ClpConstraint()
{
    type_ = 3;
    column_ = NULL;
    coefficient_ = NULL;
    numberCoefficients_ = 0;
    amplInfo_ = NULL;
}

//-------------------------------------------------------------------
// Useful Constructor
//-------------------------------------------------------------------
ClpConstraintAmpl::ClpConstraintAmpl (int row, void * amplInfo)
        : ClpConstraint()
{
    type_ = 3;
    rowNumber_ = row;
    amplInfo_ = amplInfo;
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo_;
#ifndef NDEBUG
    ASL_pfgh* asl = info->asl_;
#endif
    // warning: nlc is a macro that assumes we have a variable called asl
    assert (rowNumber_ < nlc);
    numberCoefficients_ = info->rowStart_[rowNumber_+1] - info->rowStart_[rowNumber_];
    column_ = CoinCopyOfArray(info->column_ + info->rowStart_[rowNumber_], numberCoefficients_);
    coefficient_ = new double [numberCoefficients_];;
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpConstraintAmpl::ClpConstraintAmpl (const ClpConstraintAmpl & rhs)
        : ClpConstraint(rhs)
{
    numberCoefficients_ = rhs.numberCoefficients_;
    column_ = CoinCopyOfArray(rhs.column_, numberCoefficients_);
    coefficient_ = CoinCopyOfArray(rhs.coefficient_, numberCoefficients_);
}


//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpConstraintAmpl::~ClpConstraintAmpl ()
{
    delete [] column_;
    delete [] coefficient_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpConstraintAmpl &
ClpConstraintAmpl::operator=(const ClpConstraintAmpl & rhs)
{
    if (this != &rhs) {
        delete [] column_;
        delete [] coefficient_;
        numberCoefficients_ = rhs.numberCoefficients_;
        column_ = CoinCopyOfArray(rhs.column_, numberCoefficients_);
        coefficient_ = CoinCopyOfArray(rhs.coefficient_, numberCoefficients_);
    }
    return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpConstraint * ClpConstraintAmpl::clone() const
{
    return new ClpConstraintAmpl(*this);
}

// Returns gradient
int
ClpConstraintAmpl::gradient(const ClpSimplex * model,
                            const double * solution,
                            double * gradient,
                            double & functionValue,
                            double & offset,
                            bool useScaling,
                            bool refresh) const
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo_;
    ASL_pfgh* asl = info->asl_;
    int numberColumns = n_var;;
    // If not done then do all
    if (!info->jacval_called_with_current_x_) {
        bool getStuff = eval_g(amplInfo_, numberColumns, solution, true, info->constraintValues_);
        assert (getStuff);
        getStuff = eval_jac_g(amplInfo_, numberColumns, solution, false, info->gradient_);
        assert (getStuff);
        info->jacval_called_with_current_x_ = getStuff;
    }
    if (refresh || !lastGradient_) {
        functionValue_ = info->constraintValues_[rowNumber_];
        offset_ = functionValue_; // sign??
        if (!lastGradient_)
            lastGradient_ = new double[numberColumns];
        CoinZeroN(lastGradient_, numberColumns);
        assert (!(model && model->rowScale() && useScaling));
        int i;
        int start = info->rowStart_[rowNumber_];
        assert (numberCoefficients_ == info->rowStart_[rowNumber_+1] - start);
        for (i = 0; i < numberCoefficients_; i++) {
            int iColumn = column_[i];
            double valueS = solution[iColumn];
            double valueG = info->gradient_[start+i];
            lastGradient_[iColumn] = valueG;
            offset_ -= valueS * valueG;
        }
    }
    functionValue = functionValue_;
    offset = offset_;
    memcpy(gradient, lastGradient_, numberColumns*sizeof(double));
    return 0;
}
// Resize constraint
void
ClpConstraintAmpl::resize(int newNumberColumns)
{
    abort();
}
// Delete columns in  constraint
void
ClpConstraintAmpl::deleteSome(int numberToDelete, const int * which)
{
    if (numberToDelete) {
        abort();
    }
}
// Scale constraint
void
ClpConstraintAmpl::reallyScale(const double * columnScale)
{
    abort();
}
/* Given a zeroed array sets nonlinear columns to 1.
   Returns number of nonlinear columns
*/
int
ClpConstraintAmpl::markNonlinear(char * which) const
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo_;
    ASL_pfgh* asl = info->asl_;
    int iColumn;
    int numberNon = 0;
    int nonLinear = CoinMax(nlvc, nlvo);
    for (iColumn = 0; iColumn < numberCoefficients_; iColumn++) {
        int jColumn = column_[iColumn];
        if (jColumn < nonLinear) {
            which[jColumn] = 1;
            numberNon++;
        }
    }
    return numberNon;
}
/* Given a zeroed array sets possible nonzero coefficients to 1.
   Returns number of nonzeros
*/
int
ClpConstraintAmpl::markNonzero(char * which) const
{
    int iColumn;
    for (iColumn = 0; iColumn < numberCoefficients_; iColumn++) {
        which[column_[iColumn]] = 1;
    }
    return numberCoefficients_;
}
// Number of coefficients
int
ClpConstraintAmpl::numberCoefficients() const
{
    return numberCoefficients_;
}
// Say we have new primal solution - so may need to recompute
void
ClpConstraintAmpl::newXValues()
{
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo_;
    info->conval_called_with_current_x_ = false;
    info->objval_called_with_current_x_ = false;
    info->jacval_called_with_current_x_ = false;
}
/* Load nonlinear part of problem from AMPL info
   Returns 0 if linear
   1 if quadratic objective
   2 if quadratic constraints
   3 if nonlinear objective
   4 if nonlinear constraints
   -1 on failure
*/
int
ClpSimplex::loadNonLinear(void * amplInfo, int & numberConstraints,
                          ClpConstraint ** & constraints)
{
    numberConstraints = 0;
    constraints = NULL;
    CbcAmplInfo * info = (CbcAmplInfo *) amplInfo;
    ASL_pfgh* asl = info->asl_;
    // For moment don't say quadratic
    int type = 0;
    if (nlo + nlc) {
        // nonlinear
        if (!nlc) {
            type = 3;
            delete objective_;
            objective_ = new ClpAmplObjective(amplInfo);
        } else {
            type = 4;
            numberConstraints = nlc;
            constraints = new ClpConstraint * [numberConstraints];
            if (nlo) {
                delete objective_;
                objective_ = new ClpAmplObjective(amplInfo);
            }
            for (int i = 0; i < numberConstraints; i++) {
                constraints[i] = new ClpConstraintAmpl(i, amplInfo);
            }
        }
    }
    return type;
}
#else
#include "ClpSimplex.hpp"
#include "ClpConstraint.hpp"
int
ClpSimplex::loadNonLinear(void * , int & ,
                          ClpConstraint ** & )
{
    abort();
    return 0;
}
#endif
