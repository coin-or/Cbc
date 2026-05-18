// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef ClpAbortHandler_H
#define ClpAbortHandler_H

#include <atomic>
#include "ClpEventHandler.hpp"

/** Thread-safe LP interrupt handler.
 *
 *  Attach to a ClpSimplex via passInEventHandler(). When the abort flag
 *  is set (from any thread), the next simplex iteration will stop the LP
 *  solve immediately with status 5 (stopped by event).
 *
 *  Usage:
 *  \code
 *    std::atomic<bool> abort{false};
 *    ClpAbortHandler handler(&abort);
 *    clpSimplex->passInEventHandler(&handler);
 *    // ... from another thread: abort.store(true);
 *    // LP solve will stop at next iteration
 *  \endcode
 */
class ClpAbortHandler : public ClpEventHandler {
public:
  ClpAbortHandler(const std::atomic<bool> *flag)
    : ClpEventHandler()
    , flag_(flag) {}

  ClpAbortHandler(const ClpAbortHandler &rhs)
    : ClpEventHandler(rhs)
    , flag_(rhs.flag_) {}

  ClpEventHandler *clone() const override { return new ClpAbortHandler(*this); }

  int event(Event whichEvent) override
  {
    if (whichEvent == endOfIteration && flag_ &&
        flag_->load(std::memory_order_relaxed))
      return 0; // stop LP immediately
    return -1;  // continue
  }

private:
  const std::atomic<bool> *flag_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
