/* $Id: CbcThread.hpp 1221 2009-09-03 15:43:48Z forrest $ */
// Copyright (C) 2009, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcThread_H
#define CbcThread_H

#include "CbcModel.hpp"
#include "CbcNode.hpp"
class OsiObject;
class OsiCuts;
#ifdef CBC_THREAD
class CbcThread;
// Use pthreads
#define CBC_PTHREAD
#ifdef CBC_PTHREAD
#include <pthread.h>
typedef struct {
    pthread_t	thr;
    long		status;
} Coin_pthread_t;
#endif
//#define THREAD_DEBUG 1
/** A class to encapsulate specific thread stuff
    To use another api with same style - you just need to implement
    these methods.

    At present just pthreads
 */


class CbcSpecificThread {
public:
    // Default Constructor
    CbcSpecificThread ();

    // Useful Constructor
    CbcSpecificThread (CbcSpecificThread * master, pthread_mutex_t * masterMutex);

    virtual ~CbcSpecificThread();

    // Useful stuff
    void setUsefulStuff (CbcSpecificThread * master,
                         void *& masterMutex);
    /**
       Locks a thread if parallel so that stuff like cut pool
       can be updated and/or used.
    */
    void lockThread();
    /**
       Unlocks a thread if parallel to say cut pool stuff not needed
    */
    void unlockThread();
    ///  Locks a thread for testing whether to start etc
    void lockThread2(bool doAnyway = false);
    ///  Unlocks a thread for testing whether to start etc
    void unlockThread2(bool doAnyway = false);
    /// Signal
    void signal();
    /// Timed wait in nanoseconds - if negative then seconds
    void timedWait(int time);
    /// Actually starts a thread
    void startThread(void * (*routine ) (void *), CbcThread * thread);
    /// Exits thread (called from master) - return code should be zero
    int exit();
    /// Exits thread
    void exitThread();
    /// Get status
    int status() const;
    /// Set status
    void setStatus(int value);
    //}


public: // private:
    CbcSpecificThread * basePointer_; // for getting main mutex and threadid of base
#ifdef CBC_PTHREAD
    pthread_mutex_t *masterMutex_; // for synchronizing
    pthread_mutex_t mutex2_; // for waking up threads
    pthread_cond_t condition2_; // for waking up thread
    Coin_pthread_t threadId_;
#endif
    bool locked_; // For mutex2
};
/** A class to encapsulate thread stuff */


class CbcThread {
private:
    void gutsOfDelete();
    void gutsOfCopy(const CbcThread & rhs);

public:
    // Default Constructor
    CbcThread ();

    virtual ~CbcThread();

    /// Fills in useful stuff
    void setUsefulStuff (CbcModel * model, int deterministic,
                         CbcModel * baseModel,
                         CbcThread * master,
                         void *& masterMutex);
    /**
       Locks a thread if parallel so that stuff like cut pool
       can be updated and/or used.
    */
    void lockThread();
    /**
       Unlocks a thread if parallel to say cut pool stuff not needed
    */
    void unlockThread();

    /// Returns true if locked
    inline bool isLocked() const {
        return locked_;
    }
    /** Wait for child to have return code NOT == to currentCode
        type - 0 timed wait
               1 wait
           returns true if return code changed */
    bool wait(int type, int currentCode);
    /// Just wait for so many nanoseconds
    void waitNano(int time);
    /// Signal child to carry on
    void signal();
    /// Lock from master with mutex2 and signal before lock
    void lockFromMaster();
    /// Unlock from master with mutex2 and signal after unlock
    void unlockFromMaster();
    /// Lock from thread with mutex2 and signal before lock
    void lockFromThread();
    /// Unlock from thread with mutex2 and signal after unlock
    void unlockFromThread();
    /// Exits thread (called from master) - return code should be zero
    int exit();
    /// Exits thread
    void exitThread();
    /// Waits until returnCode_ goes to zero
    void waitThread();
    /// Get status
    inline int status() const {
        return threadStuff_.status();
    }
    /// Set status
    inline void setStatus(int value) {
        threadStuff_.setStatus( value);
    }
    /// Get return code
    inline int returnCode() const {
        return returnCode_;
    }
    /// Set return code
    inline void setReturnCode(int value) {
        returnCode_ = value;
    }
    /// Get base model
    inline CbcModel * baseModel() const {
        return baseModel_;
    }
    /// Get this model
    inline CbcModel * thisModel() const {
        return thisModel_;
    }
    /// Get node
    inline CbcNode * node() const {
        return node_;
    }
    /// Set node
    inline void setNode(CbcNode * node) {
        node_ = node;
    }
    /// Get created node
    inline CbcNode * createdNode() const {
        return createdNode_;
    }
    /// Set created node
    inline void setCreatedNode(CbcNode * node) {
        createdNode_ = node;
    }
    /// Get dantzig state
    inline int dantzigState() const {
        return dantzigState_;
    }
    /// Set dantzig state
    inline void setDantzigState(int value) {
        dantzigState_ = value;
    }
    /// Get time in thread
    inline double timeInThread() const {
        return timeInThread_;
    }
    /// Increment time in thread
    inline void incrementTimeInThread(double value) {
        timeInThread_ += value;
    }
    /// Get time waiting to start
    inline double timeWaitingToStart() const {
        return timeWaitingToStart_;
    }
    /// Increment time waiting to start
    inline void incrementTimeWaitingToStart(double value) {
        timeWaitingToStart_ += value;
    }
    /// Get time locked
    inline double timeLocked() const {
        return timeLocked_;
    }
    /// Increment time locked
    inline void incrementTimeLocked(double value) {
        timeLocked_ += value;
    }
    /// Get time waiting to lock
    inline double timeWaitingToLock() const {
        return timeWaitingToLock_;
    }
    /// Increment time waiting to lock
    inline void incrementTimeWaitingToLock(double value) {
        timeWaitingToLock_ += value;
    }
    /// Get if deterministic
    inline int deterministic() const {
        return deterministic_;
    }
    /// Get maxDeleteNode
    inline int maxDeleteNode() const {
        return maxDeleteNode_;
    }
    /// Set maxDeleteNode
    inline void setMaxDeleteNode(int value) {
        maxDeleteNode_ = value;
    }
    /// Get nDeleteNode (may be fake i.e. defaultParallelIterations_)
    inline int nDeleteNode() const {
        return nDeleteNode_;
    }
    /// Set nDeleteNode (may be fake i.e. defaultParallelIterations_)
    inline void setNDeleteNode(int value) {
        nDeleteNode_ = value;
    }
    /// Clear delNode
    inline void clearDelNode() {
        delete delNode_;
        delNode_ = NULL;
    }
    /// Set fake delNode to pass across OsiCuts
    inline void fakeDelNode(CbcNode ** delNode) {
        delNode_ = delNode;
    }
    /// Get delNode
    inline CbcNode ** delNode() const {
        return delNode_;
    }
    /// Set delNode
    inline void setDelNode(CbcNode ** delNode) {
        delNode_ = delNode;
    }
    /// Get number times locked
    inline int numberTimesLocked() const {
        return numberTimesLocked_;
    }
    /// Get number times unlocked
    inline int numberTimesUnlocked() const {
        return numberTimesUnlocked_;
    }
    /// Get number of nodes this time
    inline int nodesThisTime() const {
        return nodesThisTime_;
    }
    /// Set number of nodes this time
    inline void setNodesThisTime(int value) {
        nodesThisTime_ = value;
    }
    /// Get number of iterations this time
    inline int iterationsThisTime() const {
        return iterationsThisTime_;
    }
    /// Set number of iterations this time
    inline void setIterationsThisTime(int value) {
        iterationsThisTime_ = value;
    }
    /// Get save stuff array
    inline int * saveStuff() {
        return saveStuff_;
    }
    /// Say if locked
    inline bool locked() const {
        return locked_;
    }

public: // private:
    CbcSpecificThread threadStuff_;
    CbcModel * baseModel_;
    CbcModel * thisModel_;
    CbcNode * node_; // filled in every time
    CbcNode * createdNode_; // filled in every time on return
    CbcThread * master_; // points back to master thread
    int returnCode_; // -1 available, 0 busy, 1 finished , 2??
    double timeLocked_;
    double timeWaitingToLock_;
    double timeWaitingToStart_;
    double timeInThread_;
    double timeWhenLocked_; // time when thread got lock (in seconds)
    int numberTimesLocked_;
    int numberTimesUnlocked_;
    int numberTimesWaitingToStart_;
    int saveStuff_[2];
    int dantzigState_; // 0 unset, -1 waiting to be set, 1 set
    bool locked_;
    int nDeleteNode_;
    CbcNode ** delNode_;
    int maxDeleteNode_;
    int nodesThisTime_;
    int iterationsThisTime_;
    int deterministic_;
#ifdef THREAD_DEBUG
public:
    int threadNumber_;
    int lockCount_;
#endif
};
/** Base model */


class CbcBaseModel {
public:
    // Default Constructor
    CbcBaseModel ();

    /** Constructor with model
        type -1 cuts
              0 opportunistic
              1 deterministic */
    /** Constructor with model
        type -1 cuts
              0 opportunistic
              1 deterministic */
    CbcBaseModel (CbcModel & model, int type);

    virtual ~CbcBaseModel();

    /** Stop all threads
        -1 just check all in good state
        0 actually stop
    */
    void stopThreads(int type);

    /** Wait for threads in tree
        type 0 - tree looks empty - see if any nodes outstanding
             1 - tree not empty
         2 - finish and do statistics
        returns non-zero if keep going
    */
    int waitForThreadsInTree(int type);

    /** Wait for threads n parallel cuts
        type 0 - parallel cuts
         1 - finishing parallel cuts
    */
    void waitForThreadsInCuts(int type, OsiCuts * eachCuts, int whichGenerator);

    /// Split model and do work in deterministic parallel
    void  deterministicParallel();
    /**
       Locks a thread if parallel so that stuff like cut pool
       can be updated and/or used.
    */
    inline void lockThread() {
        children_[numberThreads_].lockThread();
    }
    /**
       Unlocks a thread if parallel to say cut pool stuff not needed
    */
    inline void unlockThread() {
        children_[numberThreads_].unlockThread();
    }

    /// Returns true if locked
    inline bool isLocked() const {
        return children_[numberThreads_].locked();
    }

    /// Returns pointer to master thread
    CbcThread * masterThread() const;

    /// Returns pointer to a thread model
    inline CbcModel * model(int i) const {
        return threadModel_[i];
    }

    /// Sets Dantzig state in children
    void setDantzigState();

private:

    /// Number of children
    int numberThreads_;
    /// Child models (with base model at end)
    CbcThread * children_;
    /** type -1 cuts
              0 opportunistic
              1 deterministic */
    int type_;
    int * threadCount_;
    CbcModel ** threadModel_;
    int numberObjects_;
    OsiObject ** saveObjects_;
    int threadStats_[6];
    int defaultParallelIterations_;
    int defaultParallelNodes_;
};
#else
// Dummy threads
/** A class to encapsulate thread stuff */


class CbcThread {
public:
    // Default Constructor
    CbcThread () {}

    virtual ~CbcThread() {}

};
/** Base model */


class CbcBaseModel {
public:
    // Default Constructor (not declared here so that CbcThread.cpp not empty)
    CbcBaseModel ();

    virtual ~CbcBaseModel() {}

};
#endif

#endif

