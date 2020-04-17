
#include <math.h>
#include "calc_infect.h"

double calc_household_infect(double kappa, double omega, int HH_size, double alpha, int severe) {

  double betah=0.55; // Scaled from betah=0.4 in influenza pandemic with R0=1.6, COVID-19 R0=2.2 (Ferguson 2020)

  return(betah_scale*betah*kappa*(1+(double)severe*(omega-1))/(pow((double)HH_size,alpha)));
}

double calc_workplace_infect(int job_status, double kappa, double omega, int workplace_size, int severe, double Iw) {

  double betap[]={0.0, 1.1, 1.1, 1.1, 0.55, 0.1375} ; // Spread in all types of schools (preschool to college) is twice that of workplace
  double psi[]={0.0, 0.1, 0.2, 0.25, 0.5, 0.5} ; // Accounts for absenteeism based on severe illness. Ferguson Nature 2006

  return(betaw_scale*Iw*betap[job_status]*kappa*(1+(double)severe*(omega*psi[job_status]-1))/((double)workplace_size));
}

double calc_community_infect(double kappa, double omega, int severe, double d) {

  /* need to work on this.  Perhaps we take a random distance for each two people based on population density, number of people in county, county area, etc. */
  double zeta[]={0.1, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 0.50, 0.25, 0.25, 0.25} ; //   # Travel related parameter for community transmission. Ferguson Nature 2006
  double fd;
  double betac=0.103 ; // Scaled from betac=0.075 in influenza pandemic with R0=1.6, COVID-19 R0=2.2 (Ferguson 2020)

  fd=1/(1+pow((d/4), 3)); //kernel density function as parameterized for GB.
  return (betac_scale*betac*kappa*fd*(1+severe*(omega-1)));
}
