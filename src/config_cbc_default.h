
/***************************************************************************/
/*           HERE DEFINE THE PROJECT SPECIFIC PUBLIC MACROS                */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Version number of project */
#define CBC_VERSION "trunk"

/* Major Version number of project */
#define CBC_VERSION_MAJOR 9999

/* Minor Version number of project */
#define CBC_VERSION_MINOR 9999

/* Release Version number of project */
#define CBC_VERSION_RELEASE 9999

#ifndef CBCLIB_EXPORT
#ifdef _WIN32
/* assuming we link against a CoinUtils DLL */
#define CBCLIB_EXPORT __declspec(dllimport)
#else
#define CBCLIB_EXPORT
#endif
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
