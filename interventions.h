#if !defined(__INTERVENTIONS_H__)
#define __INTERVENTIONS_H__

enum {
  interNone, interSweMandate,
  interCaseHouseIsolation, interSchoolClosure,
  interNonEssentialClosure, interClosureAndDistancing,
  interVolWFH, interVolSelfIsolation,
  interVolWFHDistancing, interVolSelfIsolDistancing,
  interWFH50PlusDistancing, interGoogleMobilityPlusDistancing,
  interNR
};


/* Intervention type definitions */
void setup_intervention_types(double        *interIc,
                              const double **Iwp,
                              double        *interIh,
                              double        *complyI,
                              double        *tauI,
                              int           *personinter,
                              const double  *interIw[]);

#endif // if !defined(__INTERVENTIONS_H__)

/* __INTERVENTIONS_H__ */
