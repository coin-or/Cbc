
/***************************************************************************/
/*           HERE DEFINE THE PROJECT SPECIFIC PUBLIC MACROS                */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Version number of project */
#define CBC_VERSION "0.1.0"

/* Major Version number of project */
#define CBC_VERSION_MAJOR 0

/* Minor Version number of project */
#define CBC_VERSION_MINOR 1

/* Release Version number of project */
#define CBC_VERSION_RELEASE 0

#ifndef CBCLIB_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
#define CBCLIB_EXPORT __declspec(dllimport)
#else
#define CBCLIB_EXPORT
#endif
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
