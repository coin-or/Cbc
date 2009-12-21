/* $Id: CbcThread.hpp 1221 2009-09-03 15:43:48Z forrest $ */
// Copyright (C) 2009, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcThread_H
#define CbcThread_H

#include "CbcModel.hpp"
#include "CbcNode.hpp"
class OsiObject;
class OsiCuts;
#ifdef CBC_THREAD
// Threads
#include <pthread.h>
typedef struct {
    pthread_t	thr;
    long		status;
} Coin_pthread_t;

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
/** A class to encapsulate phread stuff */


class CbcPthread {
public:
    // Default Constructor
    CbcPthread ();

    // Constructor with base model
    CbcPthread (CbcModel & model, int deterministic, CbcModel * baseModel);

    // Copy constructor
    CbcPthread ( const CbcPthread &);

    virtual ~CbcPthread();

    /// Assignment operator
    CbcPthread & operator=(const CbcPthread& rhs);
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

public:
    CbcPthread * basePointer_; // for getting main mutex and threadid of base
    pthread_mutex_t mutex_; // for waking up threads
    pthread_cond_t condition_; // for waking up thread
    Coin_pthread_t threadId_;
    bool locked_;
};

/** A class to encapsulate thread stuff */


class CbcThread {
private:
    void gutsOfDelete();
    void gutsOfCopy(const CbcThread & rhs);

public:
    // Default Constructor
    CbcThread ();

    // Constructor with base model
    CbcThread (CbcModel & model, int deterministic, CbcModel * baseModel);

    // Copy constructor
    CbcThread ( const CbcThread &);

    virtual ~CbcThread();

    /// Assignment operator
    CbcThread & operator=(const CbcThread& rhs);
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

public:
    CbcModel * baseModel_;
    CbcModel * thisModel_;
    CbcNode * node_; // filled in every time
    CbcNode * createdNode_; // filled in every time on return
    Coin_pthread_t threadIdOfBase_;
    pthread_mutex_t * mutex_; // for locking data
    pthread_mutex_t * mutex2_; // for waking up threads
    CbcThread * master_; // points back to master thread
    pthread_cond_t * condition2_; // for waking up thread
    int returnCode_; // -1 available, 0 busy, 1 finished , 2??
    double timeLocked_;
    double timeWaitingToLock_;
    double timeWaitingToStart_;
    double timeInThread_;
    int numberTimesLocked_;
    int numberTimesUnlocked_;
    int numberTimesWaitingToStart_;
    int saveStuff_[2];
    int dantzigState_; // 0 unset, -1 waiting to be set, 1 set
    struct timespec absTime_;
    bool locked_;
    int nDeleteNode_;
    CbcNode ** delNode_;
    int maxDeleteNode_;
    int nodesThisTime_;
    int iterationsThisTime_;
    int deterministic_;
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
    // CbcBaseModel (CbcModel & model, int numberThreads,
    //	void *(* function) (void *),
    //	int type);
    /** Constructor with model
        type -1 cuts
              0 opportunistic
              1 deterministic */
    CbcBaseModel (CbcModel & model, int type);

    // Copy constructor
    CbcBaseModel ( const CbcBaseModel &);

    virtual ~CbcBaseModel();

    /// Assignment operator
    CbcBaseModel & operator=(const CbcBaseModel& rhs);

    /// Stop all threads
    void stopThreads();

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
        return children_[numberThreads_].locked_;
    }

    /// Returns pointer to master thread
    inline CbcThread * masterThread() const {
        return children_ + numberThreads_;
    }

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
    pthread_mutex_t mutex_main_;
    pthread_cond_t condition_main_;
    pthread_mutex_t condition_mutex_;
    Coin_pthread_t * threadId_;
    int * threadCount_;
    CbcModel ** threadModel_;
    pthread_mutex_t * mutex2_;
    pthread_cond_t * condition2_;
    int numberObjects_;
    OsiObject ** saveObjects_;
    int threadStats_[6];
    int defaultParallelIterations_;
    int defaultParallelNodes_;
};

/** Simple general method - just passed n bundles of data and a function */


class CbcSimpleThread {
public:
    // Default Constructor
    CbcSimpleThread ();

    // Constructor with stuff
    CbcSimpleThread (int numberThreads,
                     void *(* function) (void *),
                     int sizeOfData,
                     void * data);
    // Constructor with stuff (type 0 -> doHeurThread)
    CbcSimpleThread (int numberThreads,
                     int type,
                     int sizeOfData,
                     void * data);

    virtual ~CbcSimpleThread();

protected:

    /// Number of children
    int numberThreads_;
    /// data
    void * argBundle_;
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
