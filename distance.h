#if !defined(__DISTANCE_H__)
#define __DISTANCE_H__

#include <math.h>

#define deg2rad(deg) ((deg) * M_PI / 180)
#define rad2deg(rad) ((rad) * 180 / M_PI)

double distance(double lat1,
                double lon1,
                double lat2,
                double lon2,
                char   unit);

#endif /* __DISTANCE_H__ */
