/* $Id: CbcThread.cpp 1261 2009-10-30 12:45:20Z forrest $ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcConfig.h"

//#define THREAD_DEBUG
#include <string>
#include <cassert>
#include <cmath>
#include <cfloat>

#include "CbcEventHandler.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CbcThread.hpp"
#include "CbcTree.hpp"
#include "CbcHeuristic.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcModel.hpp"
#include "CbcFathom.hpp"
#include "CbcSimpleIntegerDynamicPseudoCost.hpp"
#include "ClpDualRowDantzig.hpp"
#include "OsiAuxInfo.hpp"

#include "CoinTime.hpp"
#ifdef CBC_THREAD
/// Thread functions
static void * doNodesThread(void * voidInfo);
static void * doCutsThread(void * voidInfo);
static void * doHeurThread(void * voidInfo);
// Default Constructor
CbcSpecificThread::CbcSpecificThread ()
        : basePointer_(NULL),
        masterMutex_(NULL),
        locked_(false)
{
#ifdef CBC_PTHREAD
    pthread_mutex_init(&mutex2_, NULL);
    pthread_cond_init(&condition2_, NULL);
    threadId_.status = 0;
#else
#endif
}
// Useful Constructor
CbcSpecificThread::CbcSpecificThread (CbcSpecificThread * master, pthread_mutex_t * masterMutex)
        : basePointer_(master),
        masterMutex_(masterMutex),
        locked_(false)

{
#ifdef CBC_PTHREAD
    pthread_mutex_init(&mutex2_, NULL);
    pthread_cond_init(&condition2_, NULL);
    threadId_.status = 0;
#else
#endif

}
// Useful stuff
void
CbcSpecificThread::setUsefulStuff (CbcSpecificThread * master, void *& masterMutex)

{
#ifdef CBC_PTHREAD
    basePointer_ = master;
    if (masterMutex) {
        masterMutex_ = reinterpret_cast<pthread_mutex_t *>(masterMutex);
    } else {
        // create master mutex
        masterMutex_ = new pthread_mutex_t;
        pthread_mutex_init(masterMutex_, NULL);
        masterMutex = reinterpret_cast<void *>(masterMutex_);
    }
#else
#endif
}

CbcSpecificThread::~CbcSpecificThread()
{
#ifdef CBC_PTHREAD
    pthread_mutex_destroy (&mutex2_);
    if (basePointer_ == this) {
        pthread_mutex_destroy (masterMutex_);
        delete masterMutex_;
    }
#else
#endif
}
/*
  Locks a thread if parallel so that stuff like cut pool
  can be updated and/or used.
*/
void
CbcSpecificThread::lockThread()
{
#ifdef CBC_PTHREAD
    // Use master mutex
    assert (basePointer_->masterMutex_ == masterMutex_);
    pthread_mutex_lock (masterMutex_);
#else
#endif
}
/*
  Unlocks a thread if parallel to say cut pool stuff not needed
*/
void
CbcSpecificThread::unlockThread()
{
#ifdef CBC_PTHREAD
    // Use master mutex
    pthread_mutex_unlock (masterMutex_);
#else
#endif
}
//  Locks a thread for testing whether to start etc
void
CbcSpecificThread::lockThread2(bool doAnyway)
{
    if (!locked_ || doAnyway) {
#ifdef CBC_PTHREAD
        pthread_mutex_lock (&mutex2_);
#else
#endif
        locked_ = true;
    }
}
//  Unlocks a thread for testing whether to start etc
void
CbcSpecificThread::unlockThread2(bool doAnyway)
{
    if (locked_ || doAnyway) {
#ifdef CBC_PTHREAD
        pthread_mutex_unlock (&mutex2_);
#else
#endif
        locked_ = false;
    }
}
#ifdef HAVE_CLOCK_GETTIME
inline int my_gettime(struct timespec* tp)
{
    return clock_gettime(CLOCK_REALTIME, tp);
}
#else
#ifndef _MSC_VER
inline int my_gettime(struct timespec* tp)
{
    struct timeval tv;
    int ret = gettimeofday(&tv, NULL);
    tp->tv_sec = tv.tv_sec;
    tp->tv_nsec = tv.tv_usec * 1000;
    return ret;
}
#else
inline int my_gettime(struct timespec* tp)
{
    double t = CoinGetTimeOfDay();
    tp->tv_sec = (int)floor(t);
    tp->tv_nsec = (int)((tp->tv_sec - floor(t)) / 1000000.0);
    return 0;
}
#endif
#endif
// Get time
static double getTime()
{
    struct timespec absTime2;
    my_gettime(&absTime2);
    double time2 = absTime2.tv_sec + 1.0e-9 * absTime2.tv_nsec;
    return time2;
}
// Timed wait in nanoseconds - if negative then seconds
void
CbcSpecificThread::timedWait(int time)
{

#ifdef CBC_PTHREAD
    struct timespec absTime;
    my_gettime(&absTime);
    if (time > 0) {
        absTime.tv_nsec += time;
        if (absTime.tv_nsec >= 1000000000) {
            absTime.tv_nsec -= 1000000000;
            absTime.tv_sec++;
        }
    } else {
        absTime.tv_sec -= time;
    }
    pthread_cond_timedwait(&condition2_, &mutex2_, &absTime);
#else
#endif
}
// Signal
void
CbcSpecificThread::signal()
{
#ifdef CBC_PTHREAD
    pthread_cond_signal(&condition2_);
#else
#endif
}
// Actually starts a thread
void
CbcSpecificThread::startThread(void * (*routine ) (void *), CbcThread * thread)
{
#ifdef CBC_PTHREAD
    pthread_create(&(threadId_.thr), NULL, routine, thread);
    threadId_.status = 1;
#else
#endif
}
// Exits thread (from master)
int
CbcSpecificThread::exit()
{
#ifdef CBC_PTHREAD
    pthread_cond_signal(&condition2_); // unlock
    return pthread_join(threadId_.thr, NULL);
#else
#endif
}
// Exits thread
void
CbcSpecificThread::exitThread()
{
#ifdef CBC_PTHREAD
    pthread_mutex_unlock(&mutex2_);
    pthread_exit(NULL);
#else
#endif

}
// Get status
int
CbcSpecificThread::status() const
{
#ifdef CBC_PTHREAD
    return threadId_.status;
#else
#endif
}
// Set status
void
CbcSpecificThread::setStatus(int value)
{
#ifdef CBC_PTHREAD
    threadId_.status = value;
#else
#endif
}
// Parallel heuristics
void
parallelHeuristics (int numberThreads,
                    int sizeOfData,
                    void * argBundle)
{
    Coin_pthread_t * threadId = new Coin_pthread_t [numberThreads];
    char * args = reinterpret_cast<char *>(argBundle);
    for (int i = 0; i < numberThreads; i++) {
        pthread_create(&(threadId[i].thr), NULL, doHeurThread,
                       args + i*sizeOfData);
    }
    // now wait
    for (int i = 0; i < numberThreads; i++) {
        pthread_join(threadId[i].thr, NULL);
    }
    delete [] threadId;
}
// End of specific thread stuff

/// Default constructor
CbcThread::CbcThread()
        :
        baseModel_(NULL),
        thisModel_(NULL),
        node_(NULL), // filled in every time
        createdNode_(NULL), // filled in every time on return
        returnCode_(-1), // -1 available, 0 busy, 1 finished , 2??
        timeLocked_(0.0),
        timeWaitingToLock_(0.0),
        timeWaitingToStart_(0.0),
        timeInThread_(0.0),
        numberTimesLocked_(0),
        numberTimesUnlocked_(0),
        numberTimesWaitingToStart_(0),
        dantzigState_(0), // 0 unset, -1 waiting to be set, 1 set
        locked_(false),
        nDeleteNode_(0),
        delNode_(NULL),
        maxDeleteNode_(0),
        nodesThisTime_(0),
        iterationsThisTime_(0),
        deterministic_(0)
{
}
void
CbcThread::gutsOfDelete()
{
    baseModel_ = NULL;
    thisModel_ = NULL;
    node_ = NULL;
    createdNode_ = NULL;
    delNode_ = NULL;
}
// Destructor
CbcThread::~CbcThread()
{
}
// Fills in useful stuff
void
CbcThread::setUsefulStuff (CbcModel * model, int deterministic, CbcModel * baseModel,
                           CbcThread * master,
                           void *& masterMutex)
{
    baseModel_ = baseModel;
    thisModel_ = model;
    deterministic_ = deterministic;
    threadStuff_.setUsefulStuff(&master->threadStuff_, masterMutex);
    node_ = NULL;
    createdNode_ = NULL;
    master_ = master;
    returnCode_ = -1;
    timeLocked_ = 0.0;
    timeWaitingToLock_ = 0.0;
    timeWaitingToStart_ = 0.0;
    timeInThread_ = 0.0;
    numberTimesLocked_ = 0;
    numberTimesUnlocked_ = 0;
    numberTimesWaitingToStart_ = 0;
    dantzigState_ = 0; // 0 unset, -1 waiting to be set, 1 set
    locked_ = false;
    delNode_ = NULL;
    maxDeleteNode_ = 0;
    nDeleteNode_ = 0;
    nodesThisTime_ = 0;
    iterationsThisTime_ = 0;
    if (model != baseModel) {
        // thread
        thisModel_->setInfoInChild(-3, this);
        if (deterministic_ >= 0)
            thisModel_->moveToModel(baseModel, -1);
        if (deterministic == -1)
            threadStuff_.startThread( doCutsThread, this);
        else
            threadStuff_.startThread( doNodesThread, this);
    }
}
/*
  Locks a thread if parallel so that stuff like cut pool
  can be updated and/or used.
*/
void
CbcThread::lockThread()
{
    if (!locked_) {
        double time2 = getTime();
        threadStuff_.lockThread();
        locked_ = true;
        timeWhenLocked_ = getTime();
        timeWaitingToLock_ += timeWhenLocked_ - time2;;
        numberTimesLocked_++;
#ifdef THREAD_DEBUG
        lockCount_ ++;
#if THREAD_DEBUG>1
        if (threadNumber_ == -1)
            printf("locking master %d\n", lockCount_);
        else
            printf("locking thread %d %d\n", threadNumber_, lockCount_);
#endif
    } else {
        if (threadNumber_ == -1)
            printf("master already locked %d\n", lockCount_);
        else
            printf("thread already locked %d %d\n", threadNumber_, lockCount_);
#endif
    }
}
/*
  Unlocks a thread if parallel
*/
void
CbcThread::unlockThread()
{
    if (locked_) {
        locked_ = false;
        threadStuff_.unlockThread();
        double time2 = getTime();
        timeLocked_ += time2 - timeWhenLocked_;
        numberTimesUnlocked_++;
#ifdef THREAD_DEBUG
#if THREAD_DEBUG>1
        if (threadNumber_ == -1)
            printf("unlocking master %d\n", lockCount_);
        else
            printf("unlocking thread %d %d\n", threadNumber_, lockCount_);
#endif
    } else {
        if (threadNumber_ == -1)
            printf("master already unlocked %d\n", lockCount_);
        else
            printf("thread already unlocked %d %d\n", threadNumber_, lockCount_);
#endif
    }
}
/* Wait for child to have return code NOT == to currentCode
   type - 0 timed wait
   1 wait
   returns true if return code changed */
bool
CbcThread::wait(int type, int currentCode)
{
    if (!type) {
        // just timed wait
        master_->threadStuff_.lockThread2();
        master_->threadStuff_.timedWait(1000000);
        master_->threadStuff_.unlockThread2();
    } else {
        // wait until return code changes
        while (returnCode_ == currentCode) {
            threadStuff_.signal();
            master_->threadStuff_.lockThread2();
            master_->threadStuff_.timedWait(1000000);
            master_->threadStuff_.unlockThread2();
        }
    }
    return (returnCode_ != currentCode);
}
#if 0
pthread_cond_signal(&condition2_);
-
if (!locked_)
{
    pthread_mutex_lock (&mutex2_);
    locked_ = true;
}
-
pthread_cond_timedwait(&condition2_, &mutex2_, &absTime);
-
if (locked_)
{
    pthread_mutex_unlock (&mutex2_);
    locked_ = false;
}
#endif
// Waits until returnCode_ goes to zero
void
CbcThread::waitThread()
{
    double time = getTime();
    threadStuff_.lockThread2();
    while (returnCode_) {
        threadStuff_.timedWait(-10); // 10 seconds
    }
    timeWaitingToStart_ += getTime() - time;
    numberTimesWaitingToStart_++;
}
// Just wait for so many nanoseconds
void
CbcThread::waitNano(int time)
{
    threadStuff_.lockThread2();
    threadStuff_.timedWait(time);
    threadStuff_.unlockThread2();
}
// Signal child to carry on
void
CbcThread::signal()
{
    threadStuff_.signal();
}
// Lock from master with mutex2 and signal before lock
void
CbcThread::lockFromMaster()
{
    threadStuff_.signal();
    master_->threadStuff_.lockThread2(true);
}
// Unlock from master with mutex2 and signal after unlock
void
CbcThread::unlockFromMaster()
{
    master_->threadStuff_.unlockThread2(true); // unlock anyway
    threadStuff_.signal();
}
// Lock from thread with mutex2 and signal before lock
void
CbcThread::lockFromThread()
{
    master_->threadStuff_.signal();
    threadStuff_.lockThread2();
}
// Unlock from thread with mutex2 and signal after unlock
void
CbcThread::unlockFromThread()
{
    master_->threadStuff_.signal();
    threadStuff_.unlockThread2();
}
// Exits thread (from master)
int
CbcThread::exit()
{
    return threadStuff_.exit();
}
// Exits thread
void
CbcThread::exitThread()
{
    threadStuff_.exitThread();
}
// Default constructor
CbcBaseModel::CbcBaseModel()
        :
        numberThreads_(0),
        children_(NULL),
        type_(0),
        threadCount_(NULL),
        threadModel_(NULL),
        numberObjects_(0),
        saveObjects_(NULL),
        defaultParallelIterations_(400),
        defaultParallelNodes_(2)
{
}
// Constructor with model
CbcBaseModel::CbcBaseModel (CbcModel & model,  int type)
        :
        children_(NULL),
        type_(type),
        threadCount_(NULL),
        threadModel_(NULL),
        numberObjects_(0),
        saveObjects_(NULL),
        defaultParallelIterations_(400),
        defaultParallelNodes_(2)
{
    numberThreads_ = model.getNumberThreads();
    if (numberThreads_) {
        children_ = new CbcThread [numberThreads_+1];
        // Do a partial one for base model
        void * mutex_main = NULL;
        children_[numberThreads_].setUsefulStuff(&model, type_, &model,
                children_ + numberThreads_, mutex_main);
#ifdef THREAD_DEBUG
        children_[numberThreads_].threadNumber_ = -1;
        children_[numberThreads_].lockCount_ = 0;
#endif
        threadCount_ = new int [numberThreads_];
        CoinZeroN(threadCount_, numberThreads_);
        threadModel_ = new CbcModel * [numberThreads_+1];
        memset(threadStats_, 0, sizeof(threadStats_));
        if (type_ > 0) {
            // May need for deterministic
            numberObjects_ = model.numberObjects();
            saveObjects_ = new OsiObject * [numberObjects_];
            for (int i = 0; i < numberObjects_; i++) {
                saveObjects_[i] = model.object(i)->clone();
            }
        }
        // we don't want a strategy object
        CbcStrategy * saveStrategy = model.strategy();
        model.setStrategy(NULL);
        for (int i = 0; i < numberThreads_; i++) {
            threadModel_[i] = new CbcModel(model, true);
            threadModel_[i]->synchronizeHandlers(1);
#ifdef COIN_HAS_CLP
            // Solver may need to know about model
            CbcModel * thisModel = threadModel_[i];
            CbcOsiSolver * solver =
                dynamic_cast<CbcOsiSolver *>(thisModel->solver()) ;
            if (solver)
                solver->setCbcModel(thisModel);
#endif
            children_[i].setUsefulStuff(threadModel_[i], type_, &model,
                                        children_ + numberThreads_, mutex_main);
#ifdef THREAD_DEBUG
            children_[i].threadNumber_ = i;
            children_[i].lockCount_ = 0;
#endif
        }
        model.setStrategy(saveStrategy);
    }
}
// Stop threads
void
CbcBaseModel::stopThreads(int type)
{
    if (type < 0) {
        for (int i = 0; i < numberThreads_; i++) {
            assert (children_[i].returnCode() == -1);
        }
        return;
    }
    for (int i = 0; i < numberThreads_; i++) {
        children_[i].wait(1, 0);
        assert (children_[i].returnCode() == -1);
        threadModel_[i]->setInfoInChild(-2, NULL);
        children_[i].setReturnCode( 0);
        children_[i].exit();
        children_[i].setStatus( 0);
    }
    // delete models and solvers
    for (int i = 0; i < numberThreads_; i++) {
        threadModel_[i]->setInfoInChild(type_, NULL);
        delete threadModel_[i];
    }
    delete [] children_;
    delete [] threadModel_;
    for (int i = 0; i < numberObjects_; i++)
        delete saveObjects_[i];
    delete [] saveObjects_;
    children_ = NULL;
    threadModel_ = NULL;
    saveObjects_ = NULL;
    numberObjects_ = 0;
    numberThreads_ = 0;
}
// Wait for threads in tree
int
CbcBaseModel::waitForThreadsInTree(int type)
{
    CbcModel * baseModel = children_[0].baseModel();
    int anyLeft = 0;
    // May be able to combine parts later

    if (type == 0) {
        bool locked = true;
#ifdef COIN_DEVELOP
        printf("empty\n");
#endif
        // may still be outstanding nodes
        while (true) {
            int iThread;
            for (iThread = 0; iThread < numberThreads_; iThread++) {
                if (children_[iThread].status()) {
                    if (children_[iThread].returnCode() == 0)
                        break;
                }
            }
            if (iThread < numberThreads_) {
#ifdef COIN_DEVELOP
                printf("waiting for thread %d code 0\n", iThread);
#endif
                unlockThread();
                locked = false;
                children_[iThread].wait(1, 0);
                assert(children_[iThread].returnCode() == 1);
                threadModel_[iThread]->moveToModel(baseModel, 1);
#ifdef THREAD_PRINT
                printf("off thread2 %d node %x\n", iThread, children_[iThread].node());
#endif
                children_[iThread].setNode(NULL);
                anyLeft = 1;
                assert (children_[iThread].returnCode() == 1);
                if (children_[iThread].dantzigState() == -1) {
                    // 0 unset, -1 waiting to be set, 1 set
                    children_[iThread].setDantzigState(1);
                    CbcModel * model = children_[iThread].thisModel();
                    OsiClpSolverInterface * clpSolver2
                    = dynamic_cast<OsiClpSolverInterface *> (model->solver());
                    assert (clpSolver2);
                    ClpSimplex * simplex2 = clpSolver2->getModelPtr();
                    ClpDualRowDantzig dantzig;
                    simplex2->setDualRowPivotAlgorithm(dantzig);
                }
                // say available
                children_[iThread].setReturnCode( -1);
                threadStats_[4]++;
#ifdef COIN_DEVELOP
                printf("thread %d code now -1\n", iThread);
#endif
                break;
            } else {
#ifdef COIN_DEVELOP
                printf("no threads at code 0 \n");
#endif
                // now check if any have just finished
                for (iThread = 0; iThread < numberThreads_; iThread++) {
                    if (children_[iThread].status()) {
                        if (children_[iThread].returnCode() == 1)
                            break;
                    }
                }
                if (iThread < numberThreads_) {
                    unlockThread();
                    locked = false;
                    threadModel_[iThread]->moveToModel(baseModel, 1);
#ifdef THREAD_PRINT
                    printf("off thread3 %d node %x\n", iThread, children_[iThread].node());
#endif
                    children_[iThread].setNode(NULL);
                    anyLeft = 1;
                    assert (children_[iThread].returnCode() == 1);
                    // say available
                    children_[iThread].setReturnCode( -1);
                    threadStats_[4]++;
#ifdef COIN_DEVELOP
                    printf("thread %d code now -1\n", iThread);
#endif
                    break;
                }
            }
            if (!baseModel->tree()->empty()) {
#ifdef COIN_DEVELOP
                printf("tree not empty!!!!!!\n");
#endif
                if (locked)
                    unlockThread();
                return 1;
                break;
            }
            for (iThread = 0; iThread < numberThreads_; iThread++) {
                if (children_[iThread].status()) {
                    if (children_[iThread].returnCode() != -1) {
                        printf("bad end of tree\n");
                        abort();
                    }
                }
            }
            break;
        }
#ifdef COIN_DEVELOP
        printf("finished ************\n");
#endif
        if (locked)
            unlockThread();
        return anyLeft;
    } else if (type == 1) {
        // normal
        double cutoff = baseModel->getCutoff();
        CbcNode * node = baseModel->tree()->bestNode(cutoff) ;
        // Possible one on tree worse than cutoff
        if (!node || node->objectiveValue() > cutoff)
            return 1;
        threadStats_[0]++;
        //need to think
        int iThread;
        // Start one off if any available
        for (iThread = 0; iThread < numberThreads_; iThread++) {
            if (children_[iThread].returnCode() == -1)
                break;
        }
        if (iThread < numberThreads_) {
            children_[iThread].setNode(node);
#ifdef THREAD_PRINT
            printf("empty thread %d node %x\n", iThread, children_[iThread].node());
#endif
            assert (children_[iThread].returnCode() == -1);
            // say in use
            threadModel_[iThread]->moveToModel(baseModel, 0);
            // This has to be AFTER moveToModel
            children_[iThread].setReturnCode( 0);
            children_[iThread].signal();
            threadCount_[iThread]++;
        }
        lockThread();
        // see if any finished
        for (iThread = 0; iThread < numberThreads_; iThread++) {
            if (children_[iThread].returnCode() > 0)
                break;
        }
        unlockThread();
        if (iThread < numberThreads_) {
            threadModel_[iThread]->moveToModel(baseModel, 1);
#ifdef THREAD_PRINT
            printf("off thread4 %d node %x\n", iThread, children_[iThread].node());
#endif
            children_[iThread].setNode(NULL);
            anyLeft = 1;
            assert (children_[iThread].returnCode() == 1);
            // say available
            children_[iThread].setReturnCode( -1);
            // carry on
            threadStats_[3]++;
        } else {
            // Start one off if any available
            for (iThread = 0; iThread < numberThreads_; iThread++) {
                if (children_[iThread].returnCode() == -1)
                    break;
            }
            if (iThread < numberThreads_) {
                // If any on tree get
                if (!baseModel->tree()->empty()) {
                    //node = baseModel->tree()->bestNode(cutoff) ;
                    //assert (node);
                    threadStats_[1]++;
                    return 1; // ** get another node
                }
            }
            // wait (for debug could sleep and use test)
            bool finished = false;
            while (!finished) {
                double time = getTime();
                children_[numberThreads_].wait(0, 0);
                children_[numberThreads_].incrementTimeInThread(getTime() - time);
                for (iThread = 0; iThread < numberThreads_; iThread++) {
                    if (children_[iThread].returnCode() > 0) {
                        finished = true;
                        break;
                    } else if (children_[iThread].returnCode() == 0) {
                        children_[iThread].signal(); // unlock
                    }
                }
            }
            assert (iThread < numberThreads_);
            // move information to model
            threadModel_[iThread]->moveToModel(baseModel, 1);
            anyLeft = 1;
#ifdef THREAD_PRINT
            printf("off thread %d node %x\n", iThread, children_[iThread].node());
#endif
            children_[iThread].setNode(NULL);
            assert (children_[iThread].returnCode() == 1);
            // say available
            children_[iThread].setReturnCode( -1);
        }
        // carry on
        threadStats_[2]++;
        for (int iThread = 0; iThread < numberThreads_; iThread++) {
            if (children_[iThread].status()) {
                if (children_[iThread].returnCode() != -1) {
                    anyLeft = 1;
                    break;
                }
            }
        }
        return anyLeft;
    } else if (type == 2) {
        assert (baseModel->tree()->empty());
        int i;
        for (i = 0; i < numberThreads_; i++)
            assert (children_[i].returnCode() == -1);
        // do statistics
        // Seems to be bug in CoinCpu on Linux - does threads as well despite documentation
        double time = 0.0;
        for (i = 0; i < numberThreads_; i++)
            time += children_[i].timeInThread();
        bool goodTimer = time < (baseModel->getCurrentSeconds());
        for (i = 0; i < numberThreads_; i++) {
            while (children_[i].returnCode() == 0) {
                children_[i].signal();
                double time = getTime();
                children_[numberThreads_].wait(0, 0);
                children_[numberThreads_].incrementTimeInThread(getTime() - time);
            }
            children_[i].lockFromMaster();
            threadModel_[i]->setNumberThreads(0); // say exit
            if (children_[i].deterministic() > 0)
                delete [] children_[i].delNode();
            children_[i].setReturnCode( 0);
            children_[i].unlockFromMaster();
            int returnCode;
            returnCode = children_[i].exit();
            children_[i].setStatus( 0);
            assert (!returnCode);
            //else
            threadModel_[i]->moveToModel(baseModel, 2);
            assert (children_[i].numberTimesLocked() == children_[i].numberTimesUnlocked());
            baseModel->messageHandler()->message(CBC_THREAD_STATS, baseModel->messages())
            << "Thread";
            baseModel->messageHandler()->printing(true)
            << i << threadCount_[i] << children_[i].timeWaitingToStart();
            baseModel->messageHandler()->printing(goodTimer) << children_[i].timeInThread();
            baseModel->messageHandler()->printing(false) << 0.0;
            baseModel->messageHandler()->printing(true) << children_[i].numberTimesLocked()
            << children_[i].timeLocked() << children_[i].timeWaitingToLock()
            << CoinMessageEol;
        }
        assert (children_[numberThreads_].numberTimesLocked() == children_[numberThreads_].numberTimesUnlocked());
        baseModel->messageHandler()->message(CBC_THREAD_STATS, baseModel->messages())
        << "Main thread";
        baseModel->messageHandler()->printing(false) << 0 << 0 << 0.0;
        baseModel->messageHandler()->printing(false) << 0.0;
        baseModel->messageHandler()->printing(true) << children_[numberThreads_].timeInThread();
        baseModel->messageHandler()->printing(true) << children_[numberThreads_].numberTimesLocked()
        << children_[numberThreads_].timeLocked() << children_[numberThreads_].timeWaitingToLock()
        << CoinMessageEol;
        // delete models (here in case some point to others)
        for (i = 0; i < numberThreads_; i++) {
            // make sure handler will be deleted
            threadModel_[i]->setDefaultHandler(true);
            //delete threadModel_[i];
        }
    } else {
        abort();
    }
    return 0;
}
void
CbcBaseModel::waitForThreadsInCuts(int type, OsiCuts * eachCuts,
                                   int whichGenerator)
{
    if (type == 0) {
        // cuts while doing
        bool finished = false;
        int iThread = -1;
        // see if any available
        for (iThread = 0; iThread < numberThreads_; iThread++) {
            if (children_[iThread].returnCode()) {
                finished = true;
                break;
            } else if (children_[iThread].returnCode() == 0) {
                children_[iThread].signal();
            }
        }
        while (!finished) {
            children_[numberThreads_].waitNano(1000000);
            for (iThread = 0; iThread < numberThreads_; iThread++) {
                if (children_[iThread].returnCode() > 0) {
                    finished = true;
                    break;
                } else if (children_[iThread].returnCode() == 0) {
                    children_[iThread].signal();
                }
            }
        }
        assert (iThread < numberThreads_);
        assert (children_[iThread].returnCode());
        // Use dantzigState to signal which generator
        children_[iThread].setDantzigState(whichGenerator);
        // and delNode for eachCuts
        children_[iThread].fakeDelNode(reinterpret_cast<CbcNode **> (eachCuts));
        // allow to start
        children_[iThread].setReturnCode( 0);
        children_[iThread].signal();
    } else if (type == 1) {
        // cuts - finish up
        for (int iThread = 0; iThread < numberThreads_; iThread++) {
            if (children_[iThread].returnCode() == 0) {
                bool finished = false;
                while (!finished) {
                    children_[numberThreads_].wait(0, 0);
                    if (children_[iThread].returnCode() > 0) {
                        finished = true;
                        break;
                        //#ifndef NEW_STYLE_PTHREAD
                        //} else if (children_[iThread].returnCode_ == 0) {
                        //pthread_cond_signal(&children_[iThread].threadStuff_.condition2_); // unlock
                        //#endif
                    }
                }
            }
            assert (children_[iThread].returnCode());
            // say available
            children_[iThread].setReturnCode( -1);
            //delete threadModel_[iThread]->solver();
            //threadModel_[iThread]->setSolver(NULL);
        }
    } else {
        abort();
    }
}
// Returns pointer to master thread
CbcThread *
CbcBaseModel::masterThread() const
{
    return children_ + numberThreads_;
}

// Split model and do work in deterministic parallel
void
CbcBaseModel::deterministicParallel()
{
    CbcModel * baseModel = children_[0].baseModel();
    for (int i = 0; i < numberThreads_; i++)
        threadCount_[i]++;
    int saveTreeSize = baseModel->tree()->size();
    // For now create threadModel - later modify splitModel
    CbcModel ** threadModel = new CbcModel * [numberThreads_];
    int iThread;
    for (iThread = 0; iThread < numberThreads_; iThread++)
        threadModel[iThread] = children_[iThread].thisModel();

    int nAffected = baseModel->splitModel(numberThreads_, threadModel, defaultParallelNodes_);
    // do all until finished
    for (iThread = 0; iThread < numberThreads_; iThread++) {
        // obviously tune
        children_[iThread].setNDeleteNode(defaultParallelIterations_);
    }
    // Save current state
    int iObject;
    OsiObject ** object = baseModel->objects();
    for (iObject = 0; iObject < numberObjects_; iObject++) {
        saveObjects_[iObject]->updateBefore(object[iObject]);
    }
    for (iThread = 0; iThread < numberThreads_; iThread++) {
        children_[iThread].setReturnCode( 0);
        children_[iThread].signal();
    }
    // wait
    bool finished = false;
    double time = getTime();
    while (!finished) {
        children_[numberThreads_].waitNano( 1000000); // millisecond
        finished = true;
        for (iThread = 0; iThread < numberThreads_; iThread++) {
            if (children_[iThread].returnCode() <= 0) {
                finished = false;
            }
        }
    }
    for (iThread = 0; iThread < numberThreads_; iThread++)
        children_[iThread].setReturnCode(-1);
    children_[numberThreads_].incrementTimeInThread(getTime() - time);
    // Unmark marked
    for (int i = 0; i < nAffected; i++) {
        baseModel->walkback()[i]->unmark();
    }
    int iModel;
    double scaleFactor = 1.0;
    for (iModel = 0; iModel < numberThreads_; iModel++) {
        //printf("model %d tree size %d\n",iModel,threadModel[iModel]->baseModel->tree()->size());
        if (saveTreeSize > 4*numberThreads_*defaultParallelNodes_) {
            if (!threadModel[iModel]->tree()->size()) {
                scaleFactor *= 1.05;
            }
        }
        threadModel[iModel]->moveToModel(baseModel, 11);
        // Update base model
        OsiObject ** threadObject = threadModel[iModel]->objects();
        for (iObject = 0; iObject < numberObjects_; iObject++) {
            object[iObject]->updateAfter(threadObject[iObject], saveObjects_[iObject]);
        }
    }
    if (scaleFactor != 1.0) {
        int newNumber = static_cast<int> (defaultParallelNodes_ * scaleFactor + 0.5001);
        if (newNumber*2 < defaultParallelIterations_) {
            if (defaultParallelNodes_ == 1)
                newNumber = 2;
            if (newNumber != defaultParallelNodes_) {
                char general[200];
                sprintf(general, "Changing tree size from %d to %d",
                        defaultParallelNodes_, newNumber);
                baseModel->messageHandler()->message(CBC_GENERAL,
                                                     baseModel->messages())
                << general << CoinMessageEol ;
                defaultParallelNodes_ = newNumber;
            }
        }
    }
    delete [] threadModel;
}
// Destructor
CbcBaseModel::~CbcBaseModel()
{
    delete [] threadCount_;
#if 1
    for (int i = 0; i < numberThreads_; i++)
        delete threadModel_[i];
    delete [] threadModel_;
    delete [] children_;
#endif
    for (int i = 0; i < numberObjects_; i++)
        delete saveObjects_[i];
    delete [] saveObjects_;
}
// Sets Dantzig state in children
void
CbcBaseModel::setDantzigState()
{
    for (int i = 0; i < numberThreads_; i++) {
        children_[i].setDantzigState(-1);
    }
}
static void * doNodesThread(void * voidInfo)
{
    CbcThread * stuff = reinterpret_cast<CbcThread *> (voidInfo);
    CbcModel * thisModel = stuff->thisModel();
    CbcModel * baseModel = stuff->baseModel();
    while (true) {
        stuff->waitThread();
        //printf("start node %x\n",stuff->node);
        int mode = thisModel->getNumberThreads();
        if (mode) {
            // normal
            double time2 = CoinCpuTime();
            assert (stuff->returnCode() == 0);
            if (thisModel->parallelMode() >= 0) {
                CbcNode * node = stuff->node();
                assert (node->nodeInfo());
                CbcNode * createdNode = stuff->createdNode();
                thisModel->doOneNode(baseModel, node, createdNode);
                stuff->setNode(node);
                stuff->setCreatedNode(createdNode);
                stuff->setReturnCode( 1);
            } else {
                assert (!stuff->node());
                assert (!stuff->createdNode());
                int numberIterations = stuff->nDeleteNode();
                int nDeleteNode = 0;
                int maxDeleteNode = stuff->maxDeleteNode();
                CbcNode ** delNode = stuff->delNode();
                int returnCode = 1;
                // this should be updated by heuristics strong branching etc etc
                assert (numberIterations > 0);
                thisModel->setNumberThreads(0);
                int nodesThisTime = thisModel->getNodeCount();
                int iterationsThisTime = thisModel->getIterationCount();
                int strongThisTime = thisModel->numberStrongIterations();
                thisModel->setStopNumberIterations(thisModel->getIterationCount() + numberIterations);
                int numberColumns = thisModel->getNumCols();
                int * used = CoinCopyOfArray(thisModel->usedInSolution(), numberColumns);
                int numberSolutions = thisModel->getSolutionCount();
                while (true) {
                    if (thisModel->tree()->empty()) {
                        returnCode = 1 + 1;
#ifdef CLP_INVESTIGATE_2
                        printf("%x tree empty - time %18.6f\n", thisModel, CoinGetTimeOfDay() - 1.2348e9);
#endif
                        break;
                    }
#define NODE_ITERATIONS 2
                    int nodesNow = thisModel->getNodeCount();
                    int iterationsNow = thisModel->getIterationCount();
                    int strongNow = thisModel->numberStrongIterations();
                    bool exit1 = (NODE_ITERATIONS * ((nodesNow - nodesThisTime) +
                                                     ((strongNow - strongThisTime) >> 1)) +
                                  (iterationsNow - iterationsThisTime) > numberIterations);
                    //bool exit2 =(thisModel->getIterationCount()>thisModel->getStopNumberIterations()) ;
                    //assert (exit1==exit2);
                    if (exit1 && nodesNow - nodesThisTime >= 10) {
                        // out of loop
                        //printf("out of loop\n");
#ifdef CLP_INVESTIGATE3
                        printf("%x tree %d nodes left, done %d and %d its - time %18.6f\n", thisModel,
                               thisModel->tree()->size(), nodesNow - nodesThisTime,
                               iterationsNow - iterationsThisTime, CoinGetTimeOfDay() - 1.2348e9);
#endif
                        break;
                    }
                    double cutoff = thisModel->getCutoff() ;
                    CbcNode *node = thisModel->tree()->bestNode(cutoff) ;
                    // Possible one on tree worse than cutoff
                    if (!node)
                        continue;
                    CbcNode * createdNode = NULL;
                    // Do real work of node
                    thisModel->doOneNode(NULL, node, createdNode);
                    assert (createdNode);
                    if (!createdNode->active()) {
                        delete createdNode;
                    } else {
                        // Say one more pointing to this **** postpone if marked
                        node->nodeInfo()->increment() ;
                        thisModel->tree()->push(createdNode) ;
                    }
                    if (node->active()) {
                        assert (node->nodeInfo());
                        if (node->nodeInfo()->numberBranchesLeft()) {
                            thisModel->tree()->push(node) ;
                        } else {
                            node->setActive(false);
                        }
                    } else {
                        if (node->nodeInfo()) {
                            if (!node->nodeInfo()->numberBranchesLeft())
                                node->nodeInfo()->allBranchesGone(); // can clean up
                            // So will delete underlying stuff
                            node->setActive(true);
                        }
                        if (nDeleteNode == maxDeleteNode) {
                            maxDeleteNode = (3 * maxDeleteNode) / 2 + 10;
                            stuff->setMaxDeleteNode(maxDeleteNode);
                            stuff->setDelNode(new CbcNode * [maxDeleteNode]);
                            for (int i = 0; i < nDeleteNode; i++)
                                stuff->delNode()[i] = delNode[i];
                            delete [] delNode;
                            delNode = stuff->delNode();
                        }
                        delNode[nDeleteNode++] = node;
                    }
                }
                // end of this sub-tree
                int * usedA = thisModel->usedInSolution();
                for (int i = 0; i < numberColumns; i++) {
                    usedA[i] -= used[i];
                }
                delete [] used;
                thisModel->setSolutionCount(thisModel->getSolutionCount() - numberSolutions);
                stuff->setNodesThisTime(thisModel->getNodeCount() - nodesThisTime);
                stuff->setIterationsThisTime(thisModel->getIterationCount() - iterationsThisTime);
                stuff->setNDeleteNode(nDeleteNode);
                stuff->setReturnCode( returnCode);
                thisModel->setNumberThreads(mode);
            }
            //printf("end node %x\n",stuff->node);
            stuff->unlockFromThread();
            stuff->incrementTimeInThread(CoinCpuTime() - time2);
        } else {
            // exit
            break;
        }
    }
    //printf("THREAD exiting\n");
    stuff->exitThread();
    return NULL;
}
static void * doHeurThread(void * voidInfo)
{
    typedef struct {
        double solutionValue;
        CbcModel * model;
        double * solution;
        int foundSol;
    } argBundle;
    argBundle * stuff = reinterpret_cast<argBundle *> (voidInfo);
    stuff->foundSol =
        stuff->model->heuristic(0)->solution(stuff->solutionValue,
                                             stuff->solution);
    return NULL;
}
static void * doCutsThread(void * voidInfo)
{
    CbcThread * stuff = reinterpret_cast<CbcThread *> (voidInfo);
    CbcModel * thisModel =  stuff->thisModel();
    while (true) {
        stuff->waitThread();
        //printf("start node %x\n",stuff->node);
        int mode = thisModel->getNumberThreads();
        if (mode) {
            // normal
            assert (stuff->returnCode() == 0);
            int fullScan = thisModel->getNodeCount() == 0 ? 1 : 0; //? was >0
            CbcCutGenerator * generator = thisModel->cutGenerator(stuff->dantzigState());
            generator->refreshModel(thisModel);
            OsiCuts * cuts = reinterpret_cast<OsiCuts *> (stuff->delNode());
            OsiSolverInterface * thisSolver = thisModel->solver();
            generator->generateCuts(*cuts, fullScan, thisSolver, NULL);
            stuff->setReturnCode( 1);
            stuff->unlockFromThread();
        } else {
            // exit
            break;
        }
    }
    stuff->exitThread();
    return NULL;
}
// Split up nodes - returns number of CbcNodeInfo's affected
int
CbcModel::splitModel(int numberModels, CbcModel ** model,
                     int numberNodes)
{
    int iModel;
    int i;
    for (iModel = 0; iModel < numberModels; iModel++) {
        CbcModel * otherModel = model[iModel];
        otherModel->moveToModel(this, 10);
        assert (!otherModel->tree()->size());
        otherModel->tree()->resetNodeNumbers();
        otherModel->bestPossibleObjective_ = bestPossibleObjective_;
        otherModel->sumChangeObjective1_ = sumChangeObjective1_;
        otherModel->sumChangeObjective2_ = sumChangeObjective2_;
        int numberColumns = solver_->getNumCols();
        if (otherModel->bestSolution_) {
            assert (bestSolution_);
            memcpy(otherModel->bestSolution_, bestSolution_, numberColumns*sizeof(double));
        } else if (bestSolution_) {
            otherModel->bestSolution_ = CoinCopyOfArray(bestSolution_, numberColumns);
        }
        otherModel->globalCuts_ = globalCuts_;
        otherModel->numberSolutions_ = numberSolutions_;
        otherModel->numberHeuristicSolutions_ = numberHeuristicSolutions_;
        otherModel->numberNodes_ = 1; //numberNodes_;
        otherModel->numberIterations_ = numberIterations_;
#ifdef JJF_ZERO
        if (maximumNumberCuts_ > otherModel->maximumNumberCuts_) {
            otherModel->maximumNumberCuts_ = maximumNumberCuts_;
            delete [] otherModel->addedCuts_;
            otherModel->addedCuts_ = new CbcCountRowCut * [maximumNumberCuts_];
        }
        if (maximumDepth_ > otherModel->maximumDepth_) {
            otherModel->maximumDepth_ = maximumDepth_;
            delete [] otherModel->walkback_;
            otherModel->walkback_ = new CbcNodeInfo * [maximumDepth_];
        }
#endif
        otherModel->currentNumberCuts_ = currentNumberCuts_;
        if (otherModel->usedInSolution_) {
            assert (usedInSolution_);
            memcpy(otherModel->usedInSolution_, usedInSolution_, numberColumns*sizeof(int));
        } else if (usedInSolution_) {
            otherModel->usedInSolution_ = CoinCopyOfArray(usedInSolution_, numberColumns);
        }
        /// ??? tree_;
        // Need flag (stopNumberIterations_>0?) which says don't update cut etc counts
        for (i = 0; i < numberObjects_; i++) {
            otherModel->object_[i]->updateBefore(object_[i]);
        }
        otherModel->maximumDepthActual_ = maximumDepthActual_;
        // Real cuts are in node info
        otherModel->numberOldActiveCuts_ = numberOldActiveCuts_;
        otherModel->numberNewCuts_ = numberNewCuts_;
        otherModel->numberStrongIterations_ = numberStrongIterations_;
    }
    double cutoff = getCutoff();
    int nAffected = 0;
    while (!tree_->empty()) {
        for (iModel = 0; iModel < numberModels; iModel++) {
            if (tree_->empty())
                break;
            CbcModel * otherModel = model[iModel];
            CbcNode * node = tree_->bestNode(cutoff) ;
            CbcNodeInfo * nodeInfo = node->nodeInfo();
            assert (nodeInfo);
            if (!nodeInfo->marked()) {
                //while (nodeInfo&&!nodeInfo->marked()) {
                if (nAffected == maximumDepth_) {
                    redoWalkBack();
                }
                nodeInfo->mark();
                //nodeInfo->incrementCuts(1000000);
                walkback_[nAffected++] = nodeInfo;
                //nodeInfo = nodeInfo->parent() ;
                //}
            }
            // Make node join otherModel
            OsiBranchingObject * bobj = node->modifiableBranchingObject();
            CbcBranchingObject * cbcobj = dynamic_cast<CbcBranchingObject *> (bobj);
            //assert (cbcobj);
            if (cbcobj) {
                CbcObject * object = cbcobj->object();
                assert (object);
                int position = object->position();
                assert (position >= 0);
                assert (object_[position] == object);
                CbcObject * objectNew =
                    dynamic_cast<CbcObject *> (otherModel->object_[position]);
                cbcobj->setOriginalObject(objectNew);
            }
            otherModel->tree_->push(node);
        }
        numberNodes--;
        if (!numberNodes)
            break;
    }
    return nAffected;
}
// Start threads
void
CbcModel::startSplitModel(int /*numberIterations*/)
{
    abort();
}
// Merge models
void
CbcModel::mergeModels(int /*numberModel*/, CbcModel ** /*model*/,
                      int /*numberNodes*/)
{
    abort();
}
/* Move/copy information from one model to another
   -1 - initial setup
   0 - from base model
   1 - to base model (and reset)
   2 - add in final statistics etc (and reset so can do clean destruction)
   10 - from base model (deterministic)
   11 - to base model (deterministic)
*/
void
CbcModel::moveToModel(CbcModel * baseModel, int mode)
{
    if (mode == 0) {
        setCutoff(baseModel->getCutoff());
        bestObjective_ = baseModel->bestObjective_;
        assert (!baseModel->globalCuts_.sizeRowCuts());
        numberSolutions_ = baseModel->numberSolutions_;
        stateOfSearch_ = baseModel->stateOfSearch_;
        numberNodes_ = baseModel->numberNodes_;
        numberIterations_ = baseModel->numberIterations_;
        numberFixedAtRoot_ = numberIterations_; // for statistics
        numberSolves_ = 0;
        phase_ = baseModel->phase_;
        assert (!nextRowCut_);
        nodeCompare_ = baseModel->nodeCompare_;
        tree_ = baseModel->tree_;
        assert (!subTreeModel_);
        //branchingMethod_ = NULL; // need something but what
        numberOldActiveCuts_ = baseModel->numberOldActiveCuts_;
        cutModifier_ = NULL;
        assert (!analyzeResults_);
        CbcThread * stuff = reinterpret_cast<CbcThread *> (masterThread_);
        assert (stuff);
        //if (stuff)
        stuff->setCreatedNode(NULL);
        // ?? searchStrategy_;
        searchStrategy_ = baseModel->searchStrategy_;
        stuff->saveStuff()[0] = searchStrategy_;
        stateOfSearch_ = baseModel->stateOfSearch_;
        stuff->saveStuff()[1] = stateOfSearch_;
        for (int iObject = 0 ; iObject < numberObjects_ ; iObject++) {
            CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
                dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object_[iObject]) ;
            if (dynamicObject) {
                CbcSimpleIntegerDynamicPseudoCost * baseObject =
                    dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(baseModel->object_[iObject]) ;
                assert (baseObject);
                dynamicObject->copySome(baseObject);
            }
        }
    } else if (mode == 1) {
        lockThread();
        CbcThread * stuff = reinterpret_cast<CbcThread *> (masterThread_);
        assert (stuff);
        //stateOfSearch_
        if (stuff->saveStuff()[0] != searchStrategy_) {
#ifdef COIN_DEVELOP
            printf("changing searchStrategy from %d to %d\n",
                   baseModel->searchStrategy_, searchStrategy_);
#endif
            baseModel->searchStrategy_ = searchStrategy_;
        }
        if (stuff->saveStuff()[1] != stateOfSearch_) {
#ifdef COIN_DEVELOP
            printf("changing stateOfSearch from %d to %d\n",
                   baseModel->stateOfSearch_, stateOfSearch_);
#endif
            baseModel->stateOfSearch_ = stateOfSearch_;
        }
        if (numberUpdateItems_) {
            for (int i = 0; i < numberUpdateItems_; i++) {
                CbcObjectUpdateData * update = updateItems_ + i;
                int objectNumber = update->objectNumber_;
                CbcObject * object = dynamic_cast<CbcObject *> (baseModel->object_[objectNumber]);
                if (object)
                    object->updateInformation(*update);
            }
            numberUpdateItems_ = 0;
        }
        if (eventHappened_)
            baseModel->eventHappened_ = true;
        baseModel->numberNodes_++;
        baseModel->numberIterations_ +=
            numberIterations_ - numberFixedAtRoot_;
        baseModel->numberSolves_ += numberSolves_;
        if (stuff->node())
            baseModel->tree_->push(stuff->node());
        if (stuff->createdNode())
            baseModel->tree_->push(stuff->createdNode());
        unlockThread();
    } else if (mode == 2) {
        baseModel->sumChangeObjective1_ += sumChangeObjective1_;
        baseModel->sumChangeObjective2_ += sumChangeObjective2_;
        //baseModel->numberIterations_ += numberIterations_;
        for (int iGenerator = 0; iGenerator < numberCutGenerators_; iGenerator++) {
            CbcCutGenerator * generator = baseModel->generator_[iGenerator];
            CbcCutGenerator * generator2 = generator_[iGenerator];
            generator->incrementNumberTimesEntered(generator2->numberTimesEntered());
            generator->incrementNumberCutsInTotal(generator2->numberCutsInTotal());
            generator->incrementNumberCutsActive(generator2->numberCutsActive());
            generator->incrementTimeInCutGenerator(generator2->timeInCutGenerator());
        }
        if (parallelMode() >= 0)
            nodeCompare_ = NULL;
        baseModel->maximumDepthActual_ = CoinMax(baseModel->maximumDepthActual_, maximumDepthActual_);
        baseModel->numberDJFixed_ += numberDJFixed_;
        baseModel->numberStrongIterations_ += numberStrongIterations_;
        int i;
        for (i = 0; i < 3; i++)
            baseModel->strongInfo_[i] += strongInfo_[i];
        if (parallelMode() >= 0) {
            walkback_ = NULL;
            lastNodeInfo_ = NULL;
            lastNumberCuts_ = NULL;
            lastCut_ = NULL;
            //addedCuts_ = NULL;
            tree_ = NULL;
        }
        eventHandler_ = NULL;
        delete solverCharacteristics_;
        solverCharacteristics_ = NULL;
        bool newMethod = (baseModel->branchingMethod_ && baseModel->branchingMethod_->chooseMethod());
        if (newMethod) {
            // new method - we were using base models
            numberObjects_ = 0;
            object_ = NULL;
        }
    } else if (mode == -1) {
        delete eventHandler_;
        eventHandler_ = baseModel->eventHandler_;
        assert (!statistics_);
        assert(baseModel->solverCharacteristics_);
        solverCharacteristics_ = new OsiBabSolver (*baseModel->solverCharacteristics_);
        solverCharacteristics_->setSolver(solver_);
        setMaximumNodes(COIN_INT_MAX);
        if (parallelMode() >= 0) {
            delete [] walkback_;
            //delete [] addedCuts_;
            walkback_ = NULL;
            //addedCuts_ = NULL;
            delete [] lastNodeInfo_ ;
            lastNodeInfo_ = NULL;
            delete [] lastNumberCuts_ ;
            lastNumberCuts_ = NULL;
            delete [] lastCut_ ;
            lastCut_ = NULL;
            delete tree_;
            tree_ = NULL;
            delete nodeCompare_;
            nodeCompare_ = NULL;
        } else {
            delete tree_;
            tree_ = new CbcTree();
            tree_->setComparison(*nodeCompare_) ;
        }
        continuousSolver_ = baseModel->continuousSolver_->clone();
        bool newMethod = (baseModel->branchingMethod_ && baseModel->branchingMethod_->chooseMethod());
        if (newMethod) {
            // new method uses solver - but point to base model
            // We may update an object in wrong order - shouldn't matter?
            numberObjects_ = baseModel->numberObjects_;
            if (parallelMode() >= 0) {
                object_ = baseModel->object_;
            } else {
                printf("*****WARNING - fix testosi option\n");
                object_ = baseModel->object_;
            }
        }
        int i;
        for (i = 0; i < numberHeuristics_; i++) {
            delete heuristic_[i];
            heuristic_[i] = baseModel->heuristic_[i]->clone();
            heuristic_[i]->setModelOnly(this);
        }
        for (i = 0; i < numberCutGenerators_; i++) {
            delete generator_[i];
            generator_[i] = new CbcCutGenerator(*baseModel->generator_[i]);
            // refreshModel was overkill as thought too many rows
            generator_[i]->setModel(this);
        }
    } else if (mode == 10) {
        setCutoff(baseModel->getCutoff());
        bestObjective_ = baseModel->bestObjective_;
        assert (!baseModel->globalCuts_.sizeRowCuts());
        numberSolutions_ = baseModel->numberSolutions_;
        assert (usedInSolution_);
        assert (baseModel->usedInSolution_);
        memcpy(usedInSolution_, baseModel->usedInSolution_, solver_->getNumCols()*sizeof(int));
        stateOfSearch_ = baseModel->stateOfSearch_;
        //numberNodes_ = baseModel->numberNodes_;
        //numberIterations_ = baseModel->numberIterations_;
        //numberFixedAtRoot_ = numberIterations_; // for statistics
        phase_ = baseModel->phase_;
        assert (!nextRowCut_);
        delete nodeCompare_;
        nodeCompare_ = baseModel->nodeCompare_->clone();
        tree_->setComparison(*nodeCompare_) ;
        assert (!subTreeModel_);
        //branchingMethod_ = NULL; // need something but what
        numberOldActiveCuts_ = baseModel->numberOldActiveCuts_;
        cutModifier_ = NULL;
        assert (!analyzeResults_);
        CbcThread * stuff = reinterpret_cast<CbcThread *> (masterThread_);
        assert (stuff);
        //if (stuff)
        stuff->setCreatedNode(NULL);
        // ?? searchStrategy_;
        searchStrategy_ = baseModel->searchStrategy_;
        stuff->saveStuff()[0] = searchStrategy_;
        stateOfSearch_ = baseModel->stateOfSearch_;
        stuff->saveStuff()[1] = stateOfSearch_;
        OsiObject ** baseObject = baseModel->object_;
        for (int iObject = 0 ; iObject < numberObjects_ ; iObject++) {
            object_[iObject]->updateBefore(baseObject[iObject]);
        }
        //delete [] stuff->nodeCount;
        //stuff->nodeCount = new int [baseModel->maximumDepth_+1];
    } else if (mode == 11) {
        if (parallelMode() < 0) {
            // from deterministic
            CbcThread * stuff = reinterpret_cast<CbcThread *> (masterThread_);
            assert (stuff);
            // Move solution etc
            // might as well mark all including continuous
            int numberColumns = solver_->getNumCols();
            for (int i = 0; i < numberColumns; i++) {
                baseModel->usedInSolution_[i] += usedInSolution_[i];
                //usedInSolution_[i]=0;
            }
            baseModel->numberSolutions_ += numberSolutions_;
            if (bestObjective_ < baseModel->bestObjective_ && bestObjective_ < baseModel->getCutoff()) {
                baseModel->bestObjective_ = bestObjective_ ;
                int numberColumns = solver_->getNumCols();
                if (!baseModel->bestSolution_)
                    baseModel->bestSolution_ = new double[numberColumns];
                CoinCopyN(bestSolution_, numberColumns, baseModel->bestSolution_);
                baseModel->setCutoff(getCutoff());
            }
            //stateOfSearch_
            if (stuff->saveStuff()[0] != searchStrategy_) {
#ifdef COIN_DEVELOP
                printf("changing searchStrategy from %d to %d\n",
                       baseModel->searchStrategy_, searchStrategy_);
#endif
                baseModel->searchStrategy_ = searchStrategy_;
            }
            if (stuff->saveStuff()[1] != stateOfSearch_) {
#ifdef COIN_DEVELOP
                printf("changing stateOfSearch from %d to %d\n",
                       baseModel->stateOfSearch_, stateOfSearch_);
#endif
                baseModel->stateOfSearch_ = stateOfSearch_;
            }
            int i;
            if (eventHappened_)
                baseModel->eventHappened_ = true;
            baseModel->numberNodes_ += stuff->nodesThisTime();
            baseModel->numberIterations_ += stuff->iterationsThisTime();
            double cutoff = baseModel->getCutoff();
            while (!tree_->empty()) {
                CbcNode * node = tree_->bestNode(COIN_DBL_MAX) ;
                if (node->objectiveValue() < cutoff) {
                    assert(node->nodeInfo());
                    // Make node join correctly
                    OsiBranchingObject * bobj = node->modifiableBranchingObject();
                    CbcBranchingObject * cbcobj = dynamic_cast<CbcBranchingObject *> (bobj);
                    if (cbcobj) {
                        CbcObject * object = cbcobj->object();
                        assert (object);
                        int position = object->position();
                        assert (position >= 0);
                        assert (object_[position] == object);
                        CbcObject * objectNew =
                            dynamic_cast<CbcObject *> (baseModel->object_[position]);
                        cbcobj->setOriginalObject(objectNew);
                    }
                    baseModel->tree_->push(node);
                } else {
                    delete node;
                }
            }
            for (i = 0; i < stuff->nDeleteNode(); i++) {
                //printf("CbcNode %x stuff delete\n",stuff->delNode[i]);
                delete stuff->delNode()[i];
            }
        }
    } else {
        abort();
    }
}
// Generate one round of cuts - parallel mode
int
CbcModel::parallelCuts(CbcBaseModel * master, OsiCuts & theseCuts,
                       CbcNode * /*node*/, OsiCuts & slackCuts, int lastNumberCuts)
{
    /*
      Is it time to scan the cuts in order to remove redundant cuts? If so, set
      up to do it.
    */
    int fullScan = 0 ;
    if ((numberNodes_ % SCANCUTS) == 0 || (specialOptions_&256) != 0) {
        fullScan = 1 ;
        if (!numberNodes_ || (specialOptions_&256) != 0)
            fullScan = 2;
        specialOptions_ &= ~256; // mark as full scan done
    }
    // do cuts independently
    OsiCuts * eachCuts = new OsiCuts [numberCutGenerators_];;
    int i;
    assert (master);
    for (i = 0; i < numberThreads_; i++) {
        // set solver here after cloning
        master->model(i)->solver_ = solver_->clone();
        master->model(i)->numberNodes_ = (fullScan) ? 1 : 0;
    }
    // generate cuts
    int status = 0;
    const OsiRowCutDebugger * debugger = NULL;
    bool onOptimalPath = false;
    for (i = 0; i < numberCutGenerators_; i++) {
        bool generate = generator_[i]->normal();
        // skip if not optimal and should be (maybe a cut generator has fixed variables)
        if (generator_[i]->needsOptimalBasis() && !solver_->basisIsAvailable())
            generate = false;
        if (generator_[i]->switchedOff())
            generate = false;;
        if (generate) {
            master->waitForThreadsInCuts(0, eachCuts + i, i);
        }
    }
    // wait
    master->waitForThreadsInCuts(1, eachCuts, 0);
    // Now put together
    for (i = 0; i < numberCutGenerators_; i++) {
        // add column cuts
        int numberColumnCutsBefore = theseCuts.sizeColCuts() ;
        int numberColumnCuts = eachCuts[i].sizeColCuts();
        int numberColumnCutsAfter = numberColumnCutsBefore
                                    + numberColumnCuts;
        int j;
        for (j = 0; j < numberColumnCuts; j++) {
            theseCuts.insert(eachCuts[i].colCut(j));
        }
        int numberRowCutsBefore = theseCuts.sizeRowCuts() ;
        int numberRowCuts = eachCuts[i].sizeRowCuts();
        // insert good cuts
        if (numberRowCuts) {
            int n = numberRowCuts;
            numberRowCuts = 0;
            for (j = 0; j < n; j++) {
                const OsiRowCut * thisCut = eachCuts[i].rowCutPtr(j) ;
                if (thisCut->lb() <= 1.0e10 && thisCut->ub() >= -1.0e10) {
                    theseCuts.insert(eachCuts[i].rowCut(j));
                    numberRowCuts++;
                }
            }
            if (generator_[i]->mustCallAgain() && status >= 0)
                status = 1; // say must go round
        }
        int numberRowCutsAfter = numberRowCutsBefore
                                 + numberRowCuts;
        if (numberRowCuts) {
            // Check last cut to see if infeasible
            const OsiRowCut * thisCut = theseCuts.rowCutPtr(numberRowCutsAfter - 1) ;
            if (thisCut->lb() > thisCut->ub()) {
                status = -1; // sub-problem is infeasible
                break;
            }
        }
#ifdef CBC_DEBUG
        {
            int k ;
            for (k = numberRowCutsBefore; k < numberRowCutsAfter; k++) {
                OsiRowCut thisCut = theseCuts.rowCut(k) ;
                /* check size of elements.
                   We can allow smaller but this helps debug generators as it
                   is unsafe to have small elements */
                int n = thisCut.row().getNumElements();
                const int * column = thisCut.row().getIndices();
                const double * element = thisCut.row().getElements();
                //assert (n);
                for (int i = 0; i < n; i++) {
                    double value = element[i];
                    assert(fabs(value) > 1.0e-12 && fabs(value) < 1.0e20);
                }
            }
        }
#endif
        if ((specialOptions_&1) != 0) {
            if (onOptimalPath) {
                int k ;
                for (k = numberRowCutsBefore; k < numberRowCutsAfter; k++) {
                    OsiRowCut thisCut = theseCuts.rowCut(k) ;
                    if (debugger->invalidCut(thisCut)) {
                        solver_->getRowCutDebuggerAlways()->printOptimalSolution(*solver_);
                        solver_->writeMpsNative("badCut.mps", NULL, NULL, 2);
#ifdef NDEBUG
                        printf("Cut generator %d (%s) produced invalid cut (%dth in this go)\n",
                               i, generator_[i]->cutGeneratorName(), k - numberRowCutsBefore);
                        const double *lower = getColLower() ;
                        const double *upper = getColUpper() ;
                        int numberColumns = solver_->getNumCols();
                        for (int i = 0; i < numberColumns; i++)
                            printf("%d bounds %g,%g\n", i, lower[i], upper[i]);
                        abort();
#endif
                    }
                    assert(!debugger->invalidCut(thisCut)) ;
                }
            }
        }
        /*
          The cut generator has done its thing, and maybe it generated some
          cuts.  Do a bit of bookkeeping: load
          whichGenerator[i] with the index of the generator responsible for a cut,
          and place cuts flagged as global in the global cut pool for the model.

          lastNumberCuts is the sum of cuts added in previous iterations; it's the
          offset to the proper starting position in whichGenerator.
        */
        int numberBefore =
            numberRowCutsBefore + numberColumnCutsBefore + lastNumberCuts ;
        int numberAfter =
            numberRowCutsAfter + numberColumnCutsAfter + lastNumberCuts ;
        // possibly extend whichGenerator
        resizeWhichGenerator(numberBefore, numberAfter);

        for (j = numberRowCutsBefore; j < numberRowCutsAfter; j++) {
            whichGenerator_[numberBefore++] = i ;
            const OsiRowCut * thisCut = theseCuts.rowCutPtr(j) ;
            if (thisCut->lb() > thisCut->ub())
                status = -1; // sub-problem is infeasible
            if (thisCut->globallyValid()) {
                // add to global list
                OsiRowCut newCut(*thisCut);
                newCut.setGloballyValid(true);
                newCut.mutableRow().setTestForDuplicateIndex(false);
                globalCuts_.insert(newCut) ;
            }
        }
        for (j = numberColumnCutsBefore; j < numberColumnCutsAfter; j++) {
            //whichGenerator_[numberBefore++] = i ;
            const OsiColCut * thisCut = theseCuts.colCutPtr(j) ;
            if (thisCut->globallyValid()) {
                // add to global list
                OsiColCut newCut(*thisCut);
                newCut.setGloballyValid(true);
                globalCuts_.insert(newCut) ;
            }
        }
    }
    // Add in any violated saved cuts
    if (!theseCuts.sizeRowCuts() && !theseCuts.sizeColCuts()) {
        int numberOld = theseCuts.sizeRowCuts() + lastNumberCuts;
        int numberCuts = slackCuts.sizeRowCuts() ;
        int i;
        // possibly extend whichGenerator
        resizeWhichGenerator(numberOld, numberOld + numberCuts);
        double primalTolerance;
        solver_->getDblParam(OsiPrimalTolerance, primalTolerance) ;
        for ( i = 0; i < numberCuts; i++) {
            const OsiRowCut * thisCut = slackCuts.rowCutPtr(i) ;
            if (thisCut->violated(cbcColSolution_) > 100.0*primalTolerance) {
                if (messageHandler()->logLevel() > 2)
                    printf("Old cut added - violation %g\n",
                           thisCut->violated(cbcColSolution_)) ;
                whichGenerator_[numberOld++] = -1;
                theseCuts.insert(*thisCut) ;
            }
        }
    }
    delete [] eachCuts;
    return status;
}
/*
  Locks a thread if parallel so that stuff like cut pool
  can be updated and/or used.
*/
void
CbcModel::lockThread()
{
    if (masterThread_ && (threadMode_&1) == 0)
        masterThread_->lockThread();
}
/*
  Unlocks a thread if parallel
*/
void
CbcModel::unlockThread()
{
    if (masterThread_ && (threadMode_&1) == 0)
        masterThread_->unlockThread();
}
// Returns true if locked
bool
CbcModel::isLocked() const
{
    if (masterThread_) {
        return (masterThread_->locked());
    } else {
        return true;
    }
}
// Stop a child
void
CbcModel::setInfoInChild(int type, CbcThread * info)
{
    if (type == -3) {
        // set up
        masterThread_ = info;
    } else if (type == -2) {
        numberThreads_ = 0; // signal to stop
    } else {
        // make sure message handler will be deleted
        defaultHandler_ = true;
        ownObjects_ = false;
        delete solverCharacteristics_;
        solverCharacteristics_ = NULL;
        if (type >= 0) {
            delete [] object_; // may be able to when all over to CbcThread
            for (int i = 0; i < numberCutGenerators_; i++) {
                delete generator_[i];
                generator_[i] = NULL;
                //delete virginGenerator_[i];
                //virginGenerator_[i]=NULL;
            }
            //generator_[0] = NULL;
            //delete [] generator_;
            //generator_ = NULL;
            numberCutGenerators_ = 0;
        } else {
            for (int i = 0; i < numberCutGenerators_; i++) {
                generator_[i] = NULL;
            }
        }
        object_ = NULL;
    }
}
#else
// Dummy so file not empty
// Default constructor
CbcBaseModel::CbcBaseModel() {}
#endif

