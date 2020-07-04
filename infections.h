/* Routines to manage infections */
#if !defined(__INFECTIONS_H__)
#define __INFECTIONS_H__
void place_initial_infections(char   *initial_infect_filename,
                              char   *initial_immune_filename,
                              FILE   *stats,
                              int    *initial_infections,
	                      double *tau,
                              int    *infected,
                              int    *immune,
                              int    *severe,
                              int    *symptomatic,
                              int    *county,
                              int    *num_infect,
                              int     num_counties,
                              double  symptomatic_per,
                              int     population,
                              double  dt,
                              double *lat_locale,
                              double *lon_locale,
                              int    *num_infect_county,
                              int    *num_infect_age,
                              double *age,
                              int   **county_p,
                              int    *county_size,
                              int    *locale_HH,
                              int    *HH);
#endif /* __INFECTIONS_H__ */
