#include <stdlib.h>
#include <math.h>

#include "common.h"

/* Locale related things */

int num_locale             = 0, max_locale = 0;
struct locale *locale_list = NULL;

/* Earth flattening and radius according to GRS80 taken from
  https://en.wikipedia.org/wiki/Geodetic_Reference_System_1980 */
#define FLATTENING 0.003352810681183637418
#define Earth_Radius_GRS80 6378.137

/* Earth Mean Radius according to IUGG, taken from
  https://en.wikipedia.org/wiki/Earth_radius */
#define Earth_Radius_Mean 6371.0087714

#define x(radius, lat, lon) ((radius) * \
                             sin(deg2rad(90.0 - (lat))) * cos(deg2rad((lon))))
#define y(radius, lat, lon) ((radius) * \
                             sin(deg2rad(90.0 - (lat))) * sin(deg2rad((lon))))
#define z(radius, lat) ((radius) * cos(deg2rad(90.0 - (lat))))

void add_locale(double lat, double lon, double pop_dens)
{
  if (num_locale + 1 > max_locale) {
    max_locale += 10;
    locale_list =
      (struct locale *)realloc(locale_list, max_locale * sizeof(struct locale));
  }
  locale_list[num_locale].lat   = lat;
  locale_list[num_locale].lon   = lon;
  locale_list[num_locale].rlat  = deg2rad(90.0 - lat);
  locale_list[num_locale].rlon  = deg2rad(lon);
  locale_list[num_locale].rbeta =
    atan((1 - FLATTENING) * tan(locale_list[num_locale].rlat));
  locale_list[num_locale].x           = x(Earth_Radius_Mean, lat, lon);
  locale_list[num_locale].y           = y(Earth_Radius_Mean, lat, lon);
  locale_list[num_locale].z           = z(Earth_Radius_Mean, lat);
  locale_list[num_locale].pop_density = pop_dens;

  // num_locale++;
}

/* Distance using Lambert formula, taken from
  https://en.wikipedia.org/wiki/Geographical_distance */
double locale_distance_Lambert(struct locale l1, struct locale l2)
{
  double d, ca, sca, cca, p, q, x, y, sp, sq, cp, cq;

  if ((l1.lat == l2.lat) && (l1.lon == l2.lon)) {
    return 0.0;
  }
  ca =
    acos(sin(l1.rlat) * sin(l2.rlat) + cos(l1.rlat) * cos(l2.rlat) *
         cos(l1.rlon - l2.rlon));
  sca = sin(ca / 2);
  cca = cos(ca / 2);
  p   = (l1.rbeta + l2.rbeta) / 2;
  sp  = sin(p);
  cp  = cos(p);
  q   = (l2.rbeta - l1.rbeta) / 2;
  sq  = sin(q);
  cq  = cos(q);
  x   = (ca - sin(ca)) * (sp * sp * cq * cq / (cca * cca));
  y   = (ca + sin(ca)) * (cp * cp * sq * sq / (sca * sca));
  d   = Earth_Radius_GRS80 * (ca - FLATTENING / 2 * (x + y));
  return d;
}

/* Great Circle Distance */
double locale_distance_GCD_1(struct locale l1, struct locale l2)
{
  double d;

  if ((l1.lat == l2.lat) && (l1.lon == l2.lon)) {
    return 0.0;
  }
  d =
    (l1.x * l2.x + l1.y * l2.y + l1.z *
   l2.z) / (Earth_Radius_Mean * Earth_Radius_Mean);
  return acos(d) * Earth_Radius_Mean;
}
