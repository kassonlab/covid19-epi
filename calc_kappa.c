#include <math.h>

#include "calc_kappa.h"


double calc_kappa(double t, double tau, int symptomatic, double dt, double * kappa_vals, int hosp, int icu) {

  double kappa;
  double t1;
  int t2;

  /* Determine kappa for infected person.  This is the
   * infectiousness of the person based on time since infection
   * started.  Latency period is 4.6 days.  Infection starts at
   * 5.1 days and lasts for 6 days.  Sympotmatic people are twice
   * as likely to infect others as asymptomatic.  Kappa is a log
   * normal function with mean of -0.72 and standard deviation of
   * 1.8.  From Ferguson Nature 2005
   */
  if (t-tau <= 4.6) {
    kappa=0.;
  } else if (t-tau>11.1 && hosp==0 && icu==0) {
    kappa=0.; //# Recovered or dead
  } else {
    /* First 2 lines calculates kappa on the fly, second two get precalculated kappa from array. */
    if (full_kappa) {
      t1=(log(t-tau-4.6)+0.72)/1.8;
      kappa=exp(-0.5*pow(t1,2))/((t-tau-4.6)*1.8*sqrt(2*M_PI));
    } else {
      t2=(t-tau)/dt;
      kappa=kappa_vals[t2];
    }
  }
  if (symptomatic==0) {
    kappa=kappa*0.5;
  }
  return(kappa*R0_scale);
}
