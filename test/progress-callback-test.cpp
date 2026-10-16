// Search-progress callback (Cbc_setProgressCallback / CbcEventHandler::progress):
// events must be serialised, ordered in time, monotone in both bounds, end with
// exactly one FINAL event that matches the solver's result, and be complete with
// worker threads. Maximisation and minimisation are both covered.
// This file is licensed under the Eclipse Public License (EPL).

#include "Cbc_C_Interface.h"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

namespace {
struct Monitor {
  Cbc_Model *model = nullptr;
  double sense = 1.0; // +1 minimise, -1 maximise
  std::atomic< int > inside{ 0 };
  int events = 0, incumbents = 0, finals = 0, afterFinal = 0;
  int problems = 0;
  int stopAfter = 0; // ask Cbc to stop once this many events were seen
  double lastTime = -1.0;
  double bestPrimal = 1.0e300, bestDual = -1.0e300; // minimisation sign
  double finalPrimal = 0.0, finalDual = 0.0;
  bool finalHasIncumbent = false;
  void bad(const char *what)
  {
    if (problems++ < 5)
      printf("    problem: %s\n", what);
  }
};

int CBC_LINKAGE_CB callback(Cbc_Model *model, const Cbc_ProgressInfo *info, void *data)
{
  Monitor &m = *static_cast< Monitor * >(data);
  if (m.inside.fetch_add(1) != 0)
    m.bad("callbacks overlap");
  // widen the window in which an unserialised second call would be seen
  std::this_thread::sleep_for(std::chrono::microseconds(200));
  m.events++;
  if (model != m.model)
    m.bad("wrong Cbc_Model handle");
  if (info->struct_size != (int)sizeof(Cbc_ProgressInfo))
    m.bad("struct_size");
  if (m.finals)
    m.afterFinal++;
  if (info->seconds < m.lastTime)
    m.bad("time went backwards");
  m.lastTime = info->seconds;
  if (info->has_incumbent) {
    const double p = m.sense * info->primal;
    if (p > m.bestPrimal + 1.0e-6 * (1.0 + fabs(p)))
      m.bad("incumbent got worse");
    m.bestPrimal = std::min(m.bestPrimal, p);
  }
  if (info->has_bound) {
    const double d = m.sense * info->dual;
    if (d < m.bestDual - 1.0e-6 * (1.0 + fabs(d)))
      m.bad("dual bound got worse");
    m.bestDual = std::max(m.bestDual, d);
  }
  if (info->kind == CBC_PROGRESS_INCUMBENT) {
    m.incumbents++;
    if (!info->has_incumbent)
      m.bad("INCUMBENT event without incumbent");
  }
  if (info->gap < 0.0 || info->gap > 1.0 + 1.0e-9)
    m.bad("gap out of [0,1]");
  if (!info->source || !info->source[0])
    m.bad("empty source");
  if (info->kind == CBC_PROGRESS_FINAL) {
    m.finals++;
    m.finalPrimal = info->primal;
    m.finalDual = info->dual;
    m.finalHasIncumbent = info->has_incumbent;
  }
  m.inside.fetch_sub(1);
  return m.stopAfter && m.events >= m.stopAfter;
}

// the incumbent callback shares the lock of the progress callback
int CBC_LINKAGE_CB incumbent(void *, double, int, char **, double *, void *data)
{
  Monitor &m = *static_cast< Monitor * >(data);
  if (m.inside.fetch_add(1) != 0)
    m.bad("incumbent and progress callbacks overlap");
  std::this_thread::sleep_for(std::chrono::microseconds(200));
  m.inside.fetch_sub(1);
  return 0;
}

// multi-dimensional 0-1 knapsack, big enough to branch a little
Cbc_Model *build(double sense, int n = 36)
{
  Cbc_Model *model = Cbc_newModel();
  Cbc_setLogLevel(model, 0);
  Cbc_setObjSense(model, sense < 0 ? -1.0 : 1.0);
  const int rows = n > 40 ? 10 : 4;
  unsigned seed = 12345;
  auto next = [&seed]() { seed = seed * 1103515245u + 12345u; return (seed >> 16) & 0x7fff; };
  std::vector< std::vector< double > > a(rows, std::vector< double >(n));
  std::vector< double > rhs(rows, 0.0);
  for (int j = 0; j < n; j++) {
    double obj = 20 + next() % 60;
    double w = 0;
    for (int i = 0; i < rows; i++) {
      a[i][j] = 5 + next() % 40;
      rhs[i] += a[i][j];
      w += a[i][j];
    }
    if (n > 40) // strongly correlated: hard for branch and bound
      obj = w + 10;
    std::string name = "x" + std::to_string(j);
    Cbc_addCol(model, name.c_str(), 0.0, 1.0, sense < 0 ? obj : -obj, 1, 0, nullptr, nullptr);
  }
  std::vector< int > idx(n);
  for (int j = 0; j < n; j++)
    idx[j] = j;
  for (int i = 0; i < rows; i++) {
    std::string name = "r" + std::to_string(i);
    Cbc_addRow(model, name.c_str(), n, idx.data(), a[i].data(), 'L', rhs[i] * 0.35);
  }
  return model;
}

int run(double sense, const char *threads, int stopAfter = 0)
{
  Monitor m;
  m.sense = sense;
  m.stopAfter = stopAfter;
  m.model = build(sense, stopAfter ? 70 : 36);
  Cbc_setParameter(m.model, "threads", threads);
  Cbc_setParameter(m.model, "progressInterval", "0.001");
  Cbc_setParameter(m.model, "progressMinChange", "0");
  Cbc_setProgressCallback(m.model, callback, &m);
  Cbc_addIncCallback(m.model, incumbent, &m);
  Cbc_solve(m.model);
  const bool optimal = Cbc_isProvenOptimal(m.model);
  const double obj = Cbc_getObjValue(m.model);
  if (m.events == 0)
    m.bad("no events");
  if (stopAfter) {
    if (optimal)
      m.bad("search was not stopped by the callback");
    if (m.events > stopAfter + 2 * atoi(threads) + 3)
      m.bad("too many events after the stop request");
  }
  if (m.finals != 1)
    m.bad("expected exactly one FINAL event");
  if (m.afterFinal)
    m.bad("events after FINAL");
  if (!stopAfter && m.incumbents == 0)
    m.bad("no INCUMBENT event");
  if (!stopAfter && (!m.finalHasIncumbent || fabs(m.finalPrimal - obj) > 1.0e-6 * (1.0 + fabs(obj))))
    m.bad("FINAL primal differs from the solver's objective");
  if (!stopAfter && optimal && fabs(m.sense * m.finalDual - m.sense * obj) > 1.0e-6 * (1.0 + fabs(obj)))
    m.bad("FINAL dual differs from the optimum although proven optimal");
  const bool passed = m.problems == 0;
  printf("  %s: sense=%.0f threads=%s events=%d incumbents=%d obj=%.10g optimal=%d\n",
    passed ? "ok" : "FAIL", sense, threads, m.events, m.incumbents, obj, optimal);
  Cbc_deleteModel(m.model);
  return passed ? 0 : 1;
}
} // namespace

int main()
{
  int failed = 0;
  for (double sense : { 1.0, -1.0 })
    for (const char *threads : { "0", "1", "3" })
      failed += run(sense, threads);
  for (const char *threads : { "0", "1", "3" })
    failed += run(1.0, threads, 3);
  if (failed) {
    printf("progress callback tests FAILED (%d)\n", failed);
    return 1;
  }
  printf("progress callback tests passed\n");
  return 0;
}
