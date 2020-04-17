
#ifndef CALC_INFECT_H__
#define CALC_INFECT_H__

int extern full_fd;
double extern betac_scale, betah_scale, betaw_scale;


double calc_household_infect(double, double, int, double, int);

double calc_workplace_infect(int, double, double, int, int, double);

double calc_community_infect(double, double, int, double);

#endif
