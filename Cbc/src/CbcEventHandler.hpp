/*
  Copyright (C) 2006, International Business Machines Corporation and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/

#ifndef CbcEventHandler_H
#define CbcEventHandler_H

/*! \file CbcEventHandler.hpp
    \brief Event handling for cbc

  This file contains the declaration of CbcEventHandler, used for event
  handling in cbc.

  The central method is CbcEventHandler::event(). The default semantics of
  this call are `ask for the action to take in reponse to this event'. The
  call is made at the point in the code where the event occurs (<i>e.g.</i>,
  when a solution is found, or when a node is added to or removed from the
  search tree). The return value specifies the action to perform in response
  to the event (<i>e.g.</i>, continue, or stop).

  This is a lazy class. Initially, it knows nothing about specific events,
  and returns dfltAction_ for any event. This makes for a trivial constructor
  and fast startup. The only place where the list of known events or actions
  is hardwired is in the enum definitions for CbcEvent and CbcAction,
  respectively.

  At the first call to setAction, a map is created to hold (Event,Action)
  pairs, and this map will be consulted ever after. Events not in the map
  will still return the default value.

  For serious extensions, derive a subclass and replace event() with a
  function that suits you better.  The function has access to the CbcModel
  via a pointer held in the CbcEventHandler object, and can do as much
  thinking as it likes before returning an answer.  You can also print as
  much information as you want. The model is held as a const, however, so
  you can't alter reality.

  The design of the class deliberately matches ClpEventHandler, so that other
  solvers can participate in cbc without breaking the patterns set by
  clp-specific code.

*/

#include <cstddef>
#include <map>

/* May well already be declared, but can't hurt. */

class CbcModel;

/*
  cvs/svn: $Id$
*/

/*! \class CbcEventHandler
    \brief Base class for Cbc event handling.

  Up front: We're not talking about unanticipated events here. We're talking
  about anticipated events, in the sense that the code is going to make a call
  to event() and is prepared to obey the return value that it receives.

  The general pattern for usage is as follows:
  <ol>
    <li> Create a CbcEventHandler object. This will be initialised with a set
	 of default actions for every recognised event.

    <li> Attach the event handler to the CbcModel object.

    <li> When execution reaches the point where an event occurs, call the
	 event handler as CbcEventHandler::event(the event). The return value
	 will specify what the code should do in response to the event.
  </ol>

  The return value associated with an event can be changed at any time.
*/

class CbcEventHandler {

public:
  /*! \brief Events known to cbc */

  enum CbcEvent { /*! Processing of the current node is complete. */
    node = 200,
    /*! A tree status interval has arrived. */
    treeStatus,
    /*! A solution has been found. */
    solution,
    /*! A heuristic solution has been found. */
    heuristicSolution,
    /*! A solution will be found unless user takes action (first check). */
    beforeSolution1,
    /*! A solution will be found unless user takes action (thorough check). */
    beforeSolution2,
    /*! After failed heuristic. */
    afterHeuristic,
    /*! On entry to small branch and bound. */
    smallBranchAndBound,
    /*! After a pass of heuristic. */
    heuristicPass,
    /*! When converting constraints to cuts. */
    convertToCuts,
    /*! Having generated cuts, allows user to think. */
    generatedCuts,
    /*! End of search. */
    endSearch,
    /*! Just before starting branching i.e. after root cuts. */
    afterRootCuts
  };

  /*! \brief Action codes returned by the event handler.

        Specific values are chosen to match ClpEventHandler return codes.
    */

  enum CbcAction { /*! Continue --- no action required. */
    noAction = -1,
    /*! Stop --- abort the current run at the next opportunity. */
    stop = 0,
    /*! Restart --- restart branch-and-cut search; do not undo root node
        processing.
        */
    restart,
    /*! RestartRoot --- undo root node and start branch-and-cut afresh. */
    restartRoot,
    /*! Add special cuts. */
    addCuts,
    /*! Pretend solution never happened. */
    killSolution,
    /*! Take action on modified data. */
    takeAction

  };

  /*! \brief Data type for event/action pairs */

  typedef std::map< CbcEvent, CbcAction > eaMapPair;

  /*! \name Event Processing */
  //@{

  /*! \brief Return the action to be taken for an event.

      Return the action that should be taken in response to the event passed as
      the parameter. The default implementation simply reads a return code
      from a map.
    */
  virtual CbcAction event(CbcEvent whichEvent);

  /*! \brief Return the action to be taken for an event - and modify data.

      Return the action that should be taken in response to the event passed as
      the parameter. The default implementation simply reads a return code
      from a map.
    */
  virtual CbcAction event(CbcEvent whichEvent, void *data);

  //@}

  /*! \name Constructors and destructors */
  //@{

  /*! \brief Default constructor. */

  CbcEventHandler(CbcModel *model = 0 /* was NULL but 4.6 complains */);

  /*! \brief Copy constructor. */

  CbcEventHandler(const CbcEventHandler &orig);

  /*! \brief Assignment. */

  CbcEventHandler &operator=(const CbcEventHandler &rhs);

  /*! \brief Clone (virtual) constructor. */

  virtual CbcEventHandler *clone() const;

  /*! \brief Destructor. */

  virtual ~CbcEventHandler();

  //@}

  /*! \name Set/Get methods */
  //@{

  /*! \brief Set model. */

  inline void setModel(CbcModel *model)
  {
    model_ = model;
  }

  /*! \brief Get model. */

  inline const CbcModel *getModel() const
  {
    return model_;
  }

  /*! \brief Set the default action */

  inline void setDfltAction(CbcAction action)
  {
    dfltAction_ = action;
  }

  /*! \brief Set the action code associated with an event */

  inline void setAction(CbcEvent event, CbcAction action)
  {
    if (eaMap_ == 0) {
      eaMap_ = new eaMapPair;
    }
    (*eaMap_)[event] = action;
  }

  //@}

protected:
  /*! \name Data members

       Protected (as opposed to private) to allow access by derived classes.
    */
  //@{

  /*! \brief Pointer to associated CbcModel */

  CbcModel *model_;

  /*! \brief Default action */

  CbcAction dfltAction_;

  /*! \brief Pointer to a map that holds non-default event/action pairs */

  eaMapPair *eaMap_;

  //@}
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
