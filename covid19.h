/* Defines and structure definitions */

#if !defined(__COVID19_H__)
#define __COVID19_H__

/* Things to add further on */

/*
   struct person {
    int person_idx; // Self index
    double age; // age of person
    int HH; // which household does the person live in
    int HH_sz; // Size of household
    int city; // which city does the person live in
    int county; // which county does the person live in
    int WP; // which workplace does the person go to
    int job_status; // which job status
    int infected;
    int severe;
    int infectious;
    int symptomatic;
    double long, lat;
   };

   struct household {
    int HH_idx; // Self index
    int HH_sz; // size of household
    int city;
    int county;
    int locale; // Which locale is this household in, should perhaps be pointer
       directly to the right locale
    int *persons;
    int nr_persons;
    int *infected;
    int nr_infected;
    int *infectious;
    int nr_infectious;
    double long, lat; // Probably not needed any longer
   };

   struct city  {
    int city_idx; // Self index
    int county;
    int *persons;
    int nr_persons;
    int *infected;
    int nr_infected;
    int *infectious;
    int nr_infectious;
    int *HH;
    int nr_HH;
    int *WP;
    int nr_WP;
    double long, lat;
   };

   struct county {
    int county_idx; // Keep a self index just in case we need it
    int *persons; // who lives in this county
    int nr_persons;
    int *infected; // who are infected in this county
    int nr_infected;
    int *infectious; // who are infected in this county
    int nr_infectious;
    int *HH; // which households are in this county
    int nr_HH;
    int *WP; // which workplaces are here
    int nr_WP;
    int *cities; // which cities
    int nr_cities;
   };

   struct infected {
    int infected_idx;
    int *persons;
    int nr_persons;
    int *HH;
    int nr_HH;
    int *WP;
    int nr_WP;
    int *cities;
    int nr_cities;
    int *counties;
    int nr_counties;
   };
 */

#endif /* __COVID19_H__ */
