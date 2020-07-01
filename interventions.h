#if !defined(__INTERVENTIONS_H__)
#define __INTERVENTIONS_H__

/* Intervention type definitions */
void setup_intervention_types(double       *interIc,
                              const double *Iw,
                              double       *interIh,
                              double       *complyI,
                              double       *tauI,
                              int          *personinter,
                              const double *interIw[]);

#endif // if !defined(__INTERVENTIONS_H__)

/* __INTERVENTIONS_H__ */
