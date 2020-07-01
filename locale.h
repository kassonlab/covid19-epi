#if !defined(__LOCALE_H__)
#define __LOCALE_H__

#include "common.h"

/* Locale related things */

struct locale {
  int locale_idx; // Self index

  /* Latitude and longitude */
  double lat, lon;
  double rlat, rlon; // In radians
  double rbeta;

  /* Spherical coordinates */
  double x, y, z;
  double pop_density; // Population density of this locale

  /* Things we might have a use for */

  /*
     struct county *; // Which county is this locale in
     int nPersons;
     struct person **; // Which persons are in this locale
     int nHH; // Number of households in this locale
     struct household **; // Which households are in this locale
   */
};

extern int num_locale, max_locale;
extern struct locale *locale_list;

void   add_locale(double lat,
                  double lon,
                  double pop_dens);
double locale_distance_Lambert(struct locale l1,
                               struct locale l2);
double locale_distance_GCD_1(struct locale l1,
                             struct locale l2);

#if defined(USE_LAMBERT)
# define locale_distance locale_distance_Lambert
#else // if defined(USE_LAMBERT)
# define locale_distance locale_distance_GCD_1
#endif // if defined(USE_LAMBERT)


#endif /* __LOCALE_H__ */
