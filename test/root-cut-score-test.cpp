#include "CbcCutPoolFilter.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include <cstdio>

static void add(OsiCuts &cuts, double a, double b, double rhs, double scale,
  bool mandatory = false)
{
  int indices[] = { 0, 1 };
  double coefficients[] = { a * scale, b * scale };
  OsiRowCut cut;
  cut.setRow(b ? 2 : 1, indices, coefficients);
  cut.setLb(-COIN_DBL_MAX);
  cut.setUb(rhs * scale);
  if (mandatory)
    cut.setEffectiveness(COIN_DBL_MAX);
  cuts.insert(cut);
}

int main()
{
  double x[] = { 1, 1 };
  int failures = 0;
  for (int mode : { 1, 2, 3, 4 }) {
    for (double scale : { 0.001, 1.0, 1000.0 }) {
      OsiCuts cuts;
      add(cuts, 1, 0, 0.5, scale);
      add(cuts, 1, 0, 0.6, 1);
      CbcCutPoolFilterSettings settings;
      settings.always = true;
      settings.scoreMode = mode;
      cbcFilterGeneratedCuts(cuts, 0, x, 2, 2, "test", settings);
      failures += cuts.sizeRowCuts() != 1;
      failures += cuts.rowCut(0).ub() / cuts.rowCut(0).row().getElements()[0] > 0.51;
    }
  }
  for (int mode : { 1, 2 }) {
    OsiCuts cuts;
    add(cuts, 1, 0, 0.5, 1);
    add(cuts, 1, 1, 1, 1);
    CbcCutPoolFilterSettings settings;
    settings.always = true;
    settings.maxParallelism = 1;
    settings.scoreMode = mode;
    cbcFilterGeneratedCuts(cuts, 0, x, 2, 2, "test", settings);
    failures += cuts.sizeRowCuts() != (mode == 1 ? 1 : 2);
  }
  {
    for (int mode : { 1, 3 }) {
      OsiCuts cuts;
      add(cuts, 1, 0, 0.5, 1);
      add(cuts, 1, 1, 1, 1);
      CbcCutPoolFilterSettings settings;
      settings.always = true;
      settings.scoreMode = mode;
      settings.maxParallelism = 0.9;
      cbcFilterGeneratedCuts(cuts, 0, x, 2, 2, "test", settings);
      failures += cuts.sizeRowCuts() != (mode == 1 ? 1 : 2);
    }
  }
  {
    OsiCuts cuts;
    add(cuts, 1, 0, 0.6, 1);
    add(cuts, 1, 0, 0.5, 1);
    add(cuts, 1, 0, 0.7, 1, true);
    CbcCutPoolFilterSettings settings;
    settings.always = true;
    settings.scoreMode = 1;
    std::vector< int > survivors;
    cbcFilterGeneratedCuts(cuts, 0, x, 2, 2, "test", settings, &survivors);
    failures += survivors != std::vector< int >({ 1, 2 });
    failures += cuts.sizeRowCuts() != 2;
  }
  {
    OsiCuts cuts;
    add(cuts, 1, 0, 0.5, 1);
    add(cuts, 1, 0, 0.7, 1, true);
    OsiRowCut ranged = cuts.rowCut(0);
    ranged.setLb(0);
    cuts.insert(ranged);
    CbcCutPoolFilterSettings settings;
    settings.always = true;
    settings.scoreMode = 1;
    cbcFilterGeneratedCuts(cuts, 0, x, 2, 2, "test", settings);
    failures += cuts.sizeRowCuts() != 3;
    settings.always = false;
    add(cuts, 1, 0, 0.9, 1);
    failures += cbcFilterGeneratedCuts(cuts, 0, x, 2, 2, "test", settings);
  }
  printf("root score tests: %d failures\n", failures);
  return failures != 0;
}
