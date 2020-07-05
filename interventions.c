#define Ihosp 0.25 // Accounts for infection control at hospital.

#include "interventions.h"

void setup_intervention_types(double        *interIc,
                              const double **Iwp,
                              double        *interIh,
                              double        *complyI,
                              double        *tauI,
                              int           *personinter,
                              const double  *interIw[])
{
  /* Setup intervention types. Application of interventions occurs later. */

  /* current recommendations are intervention types 1, 3, and 8. */

  /**** Introduce Intervention types.
     These are intervention types, rather than the interventions selected via
     command-line arguments. Interventions are composed from these types.
     Interventions are broken down by age group participation when necessary.
   *****/

  /* Type 0: No interventions. */
  interIc[0] = 1.0;
  static const double interIw0[6] = { 0, 1, 1, 1, 1, Ihosp };
  interIw[0]     = interIw0;
  *Iwp            = interIw0;
  interIh[0]     = 1.00;
  complyI[0]     = 1.0;
  tauI[0]        = 0;
  personinter[0] = 1;

  /* Invervention type 1: school closures, only highschools and colleges.
     No school transmission for job_status 3. Assuming no increase in community
     transmission as students would be working online at the times they would
     be in college class.*/
  interIc[1] = 1.25;
  static const double interIw1[6] = { 0, 1, 1, 0, 1, Ihosp };
  interIw[1]     = interIw1;
  interIh[1]     = 1.50;
  complyI[1]     = 1.0;
  tauI[1]        = 0;
  personinter[1] = 0;

  /* Intervention type 2: school closures of all schools.
     No school transmission for job_status 1, 2, and 3, reduction of 5%
     in workplace interactions to account for parents becoming childcare.
     Children have 50% increase in household transmission and 25% increase
     in community transmission .*/
  interIc[2] = 1.25;
  static const double interIw2[6] = { 0, 0, 0, 0, 1.0, Ihosp };
  interIw[2]     = interIw2;
  interIh[2]     = 1.50;
  complyI[2]     = 1.0;
  tauI[2]        = 0;
  personinter[2] = 0;

  /* Intervention type 3: Case isolation within household.
     1 day after symptoms start, 90% comply, household contacts remain the same,
     25% contact with community, no contact with school or workplace. */
  interIc[3] = 0.25;
  static const double interIw3[6] = { 0, 0, 0, 0, 0, Ihosp };
  interIw[3]     = interIw3;
  interIh[3]     = 1.0;
  complyI[3]     = 0.9;
  tauI[3]        = 6.1;
  personinter[3] = 0;

  /* Intervention type 4: Case isolation of entire household if one member
     becomes sick.  Same as case isoloation of single person but now includes
     all in household.  90% of symptomatic comply and 70% household members
     comply. */
  interIc[4] = 0.25;
  static const double interIw4[6] = { 0, 0, 0, 0, 0, Ihosp };
  interIw[4]     = interIw4;
  interIh[4]     = 1.5;
  complyI[4]     = 0.7;
  tauI[4]        = 0.0;
  personinter[4] = 0;

  /* Intervention type 5: Case isolation of entire household if one member
     becomes sick.  This adds for the case of a quarantined household member
     getting ill.  tauI=0 */
  interIc[5] = 0.25;
  static const double interIw5[6] = { 0, 0, 0, 0, 0, Ihosp };
  interIw[5]     = interIw5;
  interIh[5]     = 1.0;
  complyI[5]     = 0.9;
  tauI[5]        = 0.0;
  personinter[5] = 0;

  /* Intervention type 6: social distancing with school closure.
     Community contacts decrease by 75%, household comntact increase by 25%,
     70% compliance.  essential buisnesses stay open, 75% reduction in workplace
     transmission. NOTE: similar to below except with minimized social
     interaction. */
  interIc[6] = 0.25;
  static const double interIw6[6] = { 0, 0, 0, 0, 0.25, Ihosp };
  interIw[6]     = interIw6;
  interIh[6]     = 1.50;
  complyI[6]     = 1.0;
  tauI[6]        = 0;
  personinter[6] = 1;

  /* Intervention type 7: closure of schools and non-essential businesses.
     No school transmission for job_status 1, 2, and 3, reduction of 75%
     workplace interactions.  50% increase in household transmission and 50%
     increase in community transmission .*/
  interIc[7] = 1.50;
  static const double interIw7[6] = { 0, 0, 0, 0, 0.25, Ihosp };
  interIw[7]     = interIw7;
  interIh[7]     = 1.50;
  complyI[7]     = 1.0;
  tauI[7]        = 0;
  personinter[7] = 1;

  /* Intervention type 8: Social distancing of people over 70 years old.
     Reduction of 75% workplace interactions. decrease of 75% of community
     contacts, household contacts increases 25%. 80% comply*/
  interIc[8] = 0.25;
  static const double interIw8[6] = { 0, 0, 0, 0, 0.25, Ihosp };
  interIw[8]     = interIw8;
  interIh[8]     = 1.25;
  complyI[8]     = 0.8;
  tauI[8]        = 0;
  personinter[8] = 0;

  /* Intervention type 9: Social distancing of entire population.
     Reduction of 50% workplace interactions. Decrease of 75% of community
     contacts, household contacts increase 25%. Schools remain open.
     80% comply */
  interIc[9] = 0.25;
  static const double interIw9[6] = { 0, 1, 1, 0, 0.50, Ihosp };
  interIw[9]     = interIw9;
  interIh[9]     = 1.25;
  complyI[9]     = 0.8;
  tauI[9]        = 0;
  personinter[9] = 0;

  /* Intervention type 10: Social distancing of entire population as
     calculated from Google mobility data as of 2020 April 11. Reduction of
     24% workplace interactions. Decrease of 41% of community contacts.
     Schools remain open. 100% comply*/
  interIc[10] = 0.59;
  static const double interIw10[6] = { 0, 1, 1, 0, 0.76, Ihosp };
  interIw[10]     = interIw10;
  interIh[10]     = 1.00;
  complyI[10]     = 1.0;
  tauI[10]        = 0;
  personinter[10] = 1;

  /* Intervention type 11(B): Mild distancing.
     Workplace 100%, community 75%, household 100%.  100% comply */
  interIc[11] = 0.75;
  static const double interIw11[6] = { 0, 1, 1, 0, 1, Ihosp };
  interIw[11]     = interIw11;
  interIh[11]     = 1.00;
  complyI[11]     = 1.0;
  tauI[11]        = 0;
  personinter[11] = 1; // unused?

  /* Intervention type 12(C): Work from home.
     Workplace 0%, community 75%, household 150%. 100% comply */
  interIc[12] = 0.75;
  static const double interIw12[6] = { 0, 0, 0, 0, 0, Ihosp };
  interIw[12]     = interIw12;
  interIh[12]     = 1.50;
  complyI[12]     = 1.0;
  tauI[12]        = 0;
  personinter[12] = 1; // unused?

  /* Intervention type 13(D): Self-isolate without sickness.
     Workplace 0%, community 25%, household 200%. 100% comply */
  interIc[13] = 0.25;
  static const double interIw13[6] = { 0, 0, 0, 0, 0, Ihosp };
  interIw[13]     = interIw13;
  interIh[13]     = 2.00;
  complyI[13]     = 1.0;
  tauI[13]        = 0;
  personinter[13] = 1; // unused?
}
