#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

#include "CbcSppGraspGuide.hpp"

int main()
{
  const std::vector< double > rootLp{ 0.45, 0.55, 0.0 };
  const std::vector< double > roundedAttempt{ 0.0, 1.0, 0.0 };
  const std::vector< double > rawSourcePoint{ 0.4, 0.6, 0.0 };
  const std::vector< double > lower{ 0.0, 0.0, 0.0 };
  const std::vector< double > upper{ 1.0, 1.0, 1.0 };
  const std::vector< int > integerColumns{ 0, 1 };
  double rootFractionality = 0.0;
  int rootFractionalCount = 0;

  assert(CbcSppGraspGuide::preferRawSourcePoint(rootLp, rawSourcePoint, lower, upper,
    integerColumns, integerColumns, 0.8, 2, 1.0e-6, 1.0e-7,
    rootFractionality, rootFractionalCount));
  assert(std::fabs(rootFractionality - 0.9) < 1.0e-12);
  assert(rootFractionalCount == 2);

  assert(CbcSppGraspGuide::preferRawSourcePoint(rootLp, rawSourcePoint, lower, upper,
    integerColumns, integerColumns, 0.8, 2, 1.0e-6, 1.0e-7,
    rootFractionality, rootFractionalCount));
  assert(!CbcSppGraspGuide::preferRawSourcePoint(rootLp, roundedAttempt, lower, upper,
    integerColumns, integerColumns, 0.8, 2, 1.0e-6, 1.0e-7,
    rootFractionality, rootFractionalCount));
  assert(!CbcSppGraspGuide::preferRawSourcePoint(rootLp, rawSourcePoint, lower, upper,
    integerColumns, integerColumns, 0.9, 2, 1.0e-6, 1.0e-7,
    rootFractionality, rootFractionalCount));

  std::vector< double > invalidGuide(rawSourcePoint);
  invalidGuide[1] = std::numeric_limits< double >::quiet_NaN();
  assert(!CbcSppGraspGuide::preferRawSourcePoint(rootLp, invalidGuide, lower, upper,
    integerColumns, integerColumns, 0.8, 2, 1.0e-6, 1.0e-7,
    rootFractionality, rootFractionalCount));
  assert(!CbcSppGraspGuide::preferRawSourcePoint(rootLp,
    std::vector< double >{ 0.4, 0.6 }, lower, upper,
    integerColumns, integerColumns, 0.8, 2, 1.0e-6, 1.0e-7,
    rootFractionality, rootFractionalCount));
  assert(!CbcSppGraspGuide::preferRawSourcePoint(rootLp, rawSourcePoint, lower, upper,
    integerColumns, std::vector< int >{ 0 }, 0.8, 2, 1.0e-6, 1.0e-7,
    rootFractionality, rootFractionalCount));
  invalidGuide = rawSourcePoint;
  invalidGuide[0] = 1.5;
  assert(!CbcSppGraspGuide::preferRawSourcePoint(rootLp, invalidGuide, lower, upper,
    integerColumns, integerColumns, 0.8, 2, 1.0e-6, 1.0e-7,
    rootFractionality, rootFractionalCount));
  assert(!CbcSppGraspGuide::preferRawSourcePoint(rootLp, rawSourcePoint, lower, upper,
    integerColumns, integerColumns, 0.8, 1, 1.0e-6, 1.0e-7,
    rootFractionality, rootFractionalCount));

  const double inf = std::numeric_limits< double >::infinity();
  assert(CbcSppGraspGuide::unchangedBounds(
    std::vector< double >{ -inf, 0.0 },
    std::vector< double >{ inf, 1.0 },
    std::vector< double >{ -inf, 0.0 },
    std::vector< double >{ inf, 1.0 }));
  assert(!CbcSppGraspGuide::unchangedBounds(
    std::vector< double >{ -inf }, std::vector< double >{ inf },
    std::vector< double >{ -inf }, std::vector< double >{ 1.0 }));
  assert(!CbcSppGraspGuide::unchangedBounds(
    std::vector< double >{ -inf }, std::vector< double >{ inf },
    std::vector< double >{ -1.0 }, std::vector< double >{ inf }));
  return 0;
}
