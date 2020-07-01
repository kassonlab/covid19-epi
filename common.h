/* Common include file that (almost) every source file should include */

#if !defined(__COMMON_H__)
#define __COMMON_H__

#include <math.h>

#if !defined(M_PI)
# define M_PI 3.14159265358979323846
#endif // if !defined(M_PI)

/* Include files */
#include "COV_rand.h"
#include "covid19.h"
#include "distance.h"
#include "locale.h"
#include "prop_distribution.h"

#endif /* __COMMON_H__ */
