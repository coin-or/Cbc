/* test_utils.h — shared helpers for MIPster C interface tests.
 *
 * Provides fixture_path() for locating test MPS files at runtime regardless
 * of whether the tests are run from the build tree or after installation.
 *
 * Search order for fixture_path(fname):
 *   1. $MIPSTER_FIXTURE_DIR/<fname>           (runtime env override)
 *   2. FIXTURE_DIR/fixtures/<fname>           (compile-time build-tree path)
 *   3. ./fixtures/<fname>                     (cwd fallback)
 */

#ifndef MIPSTER_TEST_UTILS_H
#define MIPSTER_TEST_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIXTURE_DIR
#  define FIXTURE_DIR "."
#endif

static char _fp_buf[4096];

/* Returns a path to the fixture file <fname> (e.g. "nursesched-sprint02.mps.gz").
 * Uses a static buffer — consume before calling again. */
static const char *fixture_path(const char *fname)
{
  const char *env = getenv("MIPSTER_FIXTURE_DIR");
  if (env && env[0])
    snprintf(_fp_buf, sizeof(_fp_buf), "%s/%s", env, fname);
  else
    snprintf(_fp_buf, sizeof(_fp_buf), FIXTURE_DIR "/fixtures/%s", fname);
  return _fp_buf;
}

#endif /* MIPSTER_TEST_UTILS_H */
