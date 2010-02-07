/* $Id: CbcChooseVariable.hpp 1286 2009-11-09 23:33:07Z EdwinStraver $ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcChooseVariable_H
#define CbcChooseVariable_H

/*  These can be alternative strategies for choosing variables
    Any descendant can be passed in by setVariableChoice
*/

// FIXME: Do we want to define UP and DOWN constants (1 and -1) for clarity?
class CbcChooseVariable {
public:
    // Default Constructor
    CbcChooseVariable () {}

    virtual ~CbcChooseVariable() {}
    /** If strong branching, then only those passed in (and movement is that length)
        .  If not strong
        branching then all passed in and ignore movement.
        Returns which one chosen (or -1 if none).  way should be +1
        if branching up, -1 if down */
    virtual int chosen (const CbcModel * model, int numberToLookAt,
                        const int * which, const double * downMovement,
                        const double * upMovement, const double * solution,
                        int & way, double & value) = 0;

};

#endif

