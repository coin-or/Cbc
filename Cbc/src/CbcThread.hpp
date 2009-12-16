/* $Id: CbcThread.hpp 1221 2009-09-03 15:43:48Z forrest $ */
// Copyright (C) 2009, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcThread_H
#define CbcThread_H

#include "CbcModel.hpp"
#include "CbcNode.hpp"
class OsiObject;

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
};

#endif
