/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class that implements a conflict-based preprocessing.
 * It tries to extend set packing constraints considering
 * the conflict graph and using a greedy strategy.
 *
 * @file CglCliqueStrengthening.hpp
 * @brief Conflict-based preprocessing
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef CGLCLIQUESTRENGTHENING_HPP
#define CGLCLIQUESTRENGTHENING_HPP

#include "CoinMessageHandler.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinCliqueSet.hpp"
#include "CglConfig.h"

class OsiSolverInterface;
class CoinConflictGraph;

enum CliqueRowStatus {
    NotDominated = 1,
    Dominated = 2
};

/**
 * Auxiliary class for storing clique constraints.
 **/
class CliqueRows {
public:
  /**
   * Default constructor
   **/
  CliqueRows(size_t linesToReserve, size_t nzsToReserve);

  /**
   * Destructor
   **/
  ~CliqueRows();

  /**
   * Add a clique constraint.
   **/
  void addRow(size_t nz, const size_t els[], size_t rowIdx, CliqueRowStatus status);

  /**
   * Return the indexes of the i-th clique constraint.
   **/
  const size_t* row(size_t idxRow) const;

  /**
   * Return the original index of the i-th clique constraint.
   **/
  size_t origIdxRow(size_t idxRow) const;

  /**
   * Return the number of elements of the i-th clique constraint.
   **/
  size_t nz(size_t idxRow) const;

  /**
   * Return the status of the i-th clique constraint.
   **/
  CliqueRowStatus status(size_t idxRow) const;

  /**
   * Set the status of the i-th clique constraint.
   **/
  void setStatus(size_t idxRow, CliqueRowStatus status) const;

  /**
   * Return the number of clique constraints.
   **/
  size_t rows() const;

private:
  /**
   * Start position of each clique constraint.
   **/
  size_t *starts_;

  /**
   * Elements of the clique constraints.
   **/
  size_t *elements_;

  /**
   * Original indexes of the clique constraints.
   **/
  size_t *rowIdx_;

  /**
   * Number of clique constraints.
   **/
  size_t nRows_;

  /**
   * Status of each clique constraint.
   **/
  CliqueRowStatus *rowStatus_;
};

/**
 * Class that implements a conflict-based preprocessing.
 * It tries to extend set packing constraints considering
 * the conflict graph and using a greedy strategy.
 **/
class CGLLIB_EXPORT CglCliqueStrengthening {
public:
  /**
   * Default constructor
   **/
  CglCliqueStrengthening(OsiSolverInterface *model, CoinMessageHandler *dhandler = NULL);

  /**
   * Destructor
   **/
  ~CglCliqueStrengthening();

  /**
   * Tries to strengthen the set packing constraints of the
   * model. After strengthening (extending), dominated
   * constraints are removed (clique merging).
   *
   * @param extMethod Extension method: 0 = no extension;1 = random;
   * 2 = max degree; 3 = max modified degree;
   * 4 = reduced cost (inversely proportional);
   * 5 = reduced cost (inversely proportional) + modified degree.
   **/
  void strengthenCliques(size_t extMethod = 4);

  /**
   * Tries to strengthen the set packing constraints of
   * model whose indexes are in array rows. After
   * strengthening (extending), dominated constraints are
   * removed (clique merging).
   *
   * @param n number of rows to be strengthened
   * @param rows rows to be strengthened
   * @param extMethod Extension method: 0 = no extension;
   * 1 = random; 2 = max degree; 3 = max modified degree;
   * 4 = reduced cost (inversely proportional);
   * 5 = reduced cost (inversely proportional) + modified degree.
   **/
  void strengthenCliques(size_t n, const size_t rows[], size_t extMethod = 4);

  /**
   * Return the number of set packing constraints extended.
   **/
  int constraintsExtended() const { return nExtended_; }

  /**
   * Return the number of set packing constraints dominated
   * by the extended constraints.
   **/
  int constraintsDominated() const { return nDominated_; }

  /**
   * Set maximum wall-clock seconds for strengthening (0 = no limit).
   **/
  void setMaximumSeconds(double seconds) { maxSeconds_ = seconds; }

  /**
   * Pass in Message handler (not deleted at end)
   **/
  void passInMessageHandler(CoinMessageHandler * handler);

  /**
   * Set language
   **/
  void newLanguage(CoinMessages::Language language);

  /**
   * New language
   **/
  inline void setLanguage(CoinMessages::Language language)
  {newLanguage(language);}

  /**
   * Return message handler
   **/
  inline CoinMessageHandler * messageHandler() const
  {return handler_;}

  /**
   * Return messages
   **/
  inline CoinMessages messages()
  {return messages_;}

  /**
   * Return a pointer to messages
   **/
  inline CoinMessages * messagesPointer()
  {return &messages_;}

private:
  /**
   * Detect clique constraints in the MILP.
   **/
  void detectCliqueRows();

  /**
   * Compute the occurrences of the variables in the clique constraints.
   **/
  void fillCliquesByColumn();

  /**
   * Try to extend the clique constraints.
   **/
  void cliqueExtension(size_t extMethod, CoinCliqueSet *newCliques);

  /**
   * Try to extend the clique constraints, considering specific
   * clique constraints.
   **/
  void cliqueExtension(size_t extMethod, CoinCliqueSet *newCliques, size_t n, const size_t rows[]);

  /**
   * Fill and return the reduced costs of the variables.
   **/
  double* getReducedCost();

  /**
   * Check if a clique constraint dominates other clique constraints
   * stored in cliqueRows_.
   **/
  void checkDominance(const size_t *extClqEl, size_t extClqSize, bool *ivRow, bool *ivCol);

  /**
   * Remove dominated constraints.
   **/
  void removeDominatedRows();

  /**
   * Add the extended clique constraints.
   **/
  void addStrongerCliques(const CoinCliqueSet *newCliques);

  /**
   * A pointer to the MILP structure.
   **/
  OsiSolverInterface *model_;

  /**
   * A pointer to the conflict graph.
   **/
  const CoinConflictGraph *cgraph_;

  /**
   * Set of clique constraints
   **/
  CliqueRows *cliqueRows_;

  /**
   * Number of clique constraints that each variable appears.
   **/
  size_t *nColClqs_;

  /**
   * Indexes of the clique constraints that each variable appears.
   **/
  size_t **colClqs_;

  /**
   * Names of the clique constraints.
   **/
  OsiSolverInterface::OsiNameVec rowClqNames_;

  /**
   * Stores the original index of each constraint in cliqueRows_.
   * If a constraint is not a clique constraint, its index is invalid.
   **/
  size_t *posInClqRows_;

  /**
   * Number of extended constraints.
   **/
  int nExtended_;

  /**
   * Number of dominated constraints.
   **/
  int nDominated_;

  /**
   * Maximum wall-clock seconds allowed for strengthening (0 = no limit).
   **/
  double maxSeconds_;

  /**
   * Wall-clock time at the start of strengthenCliques().
   **/
  double startTime_;

  /**
   * Scratch space for tracking dirty row indices in checkDominance.
   **/
  size_t *dirtyRows_;

  /**
   * Message handler
   **/
  CoinMessageHandler * handler_;

  /**
   * Flag to say if handler_ is the default handler.
   * The default handler is deleted when the model
   * is deleted. Other handlers (supplied by the client)
   * will not be deleted.
   **/
  bool defaultHandler_;

  /**
   * Messages
   **/
  CoinMessages messages_;
};


#endif //CGLCLIQUESTRENGTHENING_HPP
