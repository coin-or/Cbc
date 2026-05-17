/* CglConfig.h — static configuration for Cgl merged into Cbc.
 * All required dependencies (CoinUtils, Clp, OsiClp) are always present.
 */
#ifndef __CGLCONFIG_H__
#define __CGLCONFIG_H__

/* Symbol visibility — empty for static builds */
#define CGLLIB_EXPORT

/* Available dependencies (always true in this build) */
#define CGL_HAS_COINUTILS 1
#define CGL_HAS_CLP 1
#define CGL_HAS_OSICLP 1

/* Version */
#define CGL_VERSION "devel"
#define CGL_VERSION_MAJOR 9999
#define CGL_VERSION_MINOR 9999
#define CGL_VERSION_RELEASE 9999

#endif /* __CGLCONFIG_H__ */
