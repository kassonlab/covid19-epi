#include <cuda_runtime.h>
#include <thrust/copy.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "locale.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* Earth flattening and radius according to GRS80 taken from https://en.wikipedia.org/wiki/Geodetic_Reference_System_1980 */
#define FLATTENING 0.003352810681183637418
#define Earth_Radius_GRS80 6378.137

/* Earth Mean Radius according to IUGG, taken from https://en.wikipedia.org/wiki/Earth_radius */
#define Earth_Radius_Mean 6371.0087714

#define x(radius, lat, lon) ((radius)*sin(deg2rad(90.0-(lat)))*cos(deg2rad((lon))))
#define y(radius, lat, lon) ((radius)*sin(deg2rad(90.0-(lat)))*sin(deg2rad((lon))))
#define z(radius, lat) ((radius)*cos(deg2rad(90.0-(lat))))

#undef locale_distance
#if defined(USE_LAMBERT)
#define locale_distance d_locale_distance_Lambert
#else
#define locale_distance d_locale_distance_GCD_1
#endif

/* Distance using Lambert formula, taken from https://en.wikipedia.org/wiki/Geographical_distance */
static __device__ double d_locale_distance_Lambert(struct locale l1, struct locale l2)
{
    double d, ca, sca, cca, p, q, x, y, sp, sq, cp, cq;

    if (l1.lat == l2.lat && l1.lon == l2.lon) {
        return 0.0;
    }
    ca = acos(sin(l1.rlat) * sin(l2.rlat) + cos(l1.rlat) * cos(l2.rlat) * cos(l1.rlon - l2.rlon));
    sca = sin(ca/2);
    cca = cos(ca/2);
    p = (l1.rbeta+l2.rbeta)/2;
    sp = sin(p);
    cp = cos(p);
    q = (l2.rbeta-l1.rbeta)/2;
    sq = sin(q);
    cq = cos(q);
    x = (ca-sin(ca))*(sp*sp*cq*cq/(cca*cca));
    y = (ca+sin(ca))*(cp*cp*sq*sq/(sca*sca));
    d = Earth_Radius_GRS80*(ca - FLATTENING/2*(x+y));
    return d;
}

/* Great Circle Distance */
static __device__ double d_locale_distance_GCD_1(struct locale l1, struct locale l2)
{
    double d;
    if (l1.lat == l2.lat && l1.lon == l2.lon) {
        return 0.0;
    }
    d = (l1.x * l2.x + l1.y * l2.y + l1.z * l2.z) / (Earth_Radius_Mean*Earth_Radius_Mean);
    return acos(d)*Earth_Radius_Mean;
}

#define deg2rad(deg) (deg * M_PI / 180)
#define rad2deg(rad) (rad * 180 / M_PI)

static __device__ double calc_community_infect(double kappa, double omega, int severe, double d, double betac_scale) {

	/* need to work on this.  Perhaps we take a random distance for each two people based on population density, number of people in county, county area, etc. */
	double zeta[]={0.1, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.75, 0.50, 0.25, 0.25, 0.25} ; //   # Travel related parameter for community transmission. Ferguson Nature 2006
	double fd;
	double betac=0.103 ; // Scaled from betac=0.075 in influenza pandemic with R0=1.6, COVID-19 R0=2.2 (Ferguson 2020)

	fd=1/(1+pow((d/4), 3)); //kernel density function as parameterized for GB.
	return (betac_scale*betac*kappa*fd*(1+severe*(omega-1)));
}

static __device__ double d_distance(double lat1, double lon1, double lat2, double lon2, char unit) {
    double theta, dist;
    if ((lat1 == lat2) && (lon1 == lon2)) {
      return 0;
    }
    else {
      theta = lon1 - lon2;
  //    dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
      double ang1,ang2;
      ang1 = deg2rad(lat1);
      ang2 = deg2rad(lat2);
      dist = cos(ang1) * cos(ang2) * ( 1.0 + cos(deg2rad(theta)) ) - cos(ang1 + ang2);
      dist = acos(dist);
      dist = rad2deg(dist);
      dist = dist * 60 * 1.1515;
      switch(unit) {
        case 'M':
          break;
        case 'K':
          dist = dist * 1.609344;
          break;
        case 'N':
          dist = dist * 0.8684;
          break;
      }
      return (dist);
    }
}

static __device__ double calc_kappa(double t, double tau, int symptomatic, double dt, double const* kappa_vals, int hosp, int icu, int full_kappa, double R0_scale) {

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

struct LoopInvariantData {
    thrust::device_vector<int> infectious;
    thrust::device_vector<int> intervene;
    thrust::device_vector<double> tau;
    thrust::device_vector<double> tauI;
    thrust::device_vector<double> interIc;
    thrust::device_vector<int> symptomatic;
    thrust::device_vector<double> kappa_vals;
    thrust::device_vector<int> hosp_pop;
    thrust::device_vector<int> icu_pop;
    thrust::device_vector<double> lat_locale;
    thrust::device_vector<double> lon_locale;
    thrust::device_vector<int> locale_HH;
    thrust::device_vector<int> HH;
    thrust::device_vector<struct locale> locale_list;
    thrust::device_vector<int> severe;
    thrust::device_vector<double> kappa_cache;
};


static __global__ void locale_infectious_step_kernel(
    int j,
    int num_infectious,
    double Ic,
    double t,
    double dt,
    double omega,
    int full_kappa,
    double R0_scale,
    double betac_scale,

    int const* infectious,
    int const* intervene,
    double const* tau,
    double const* tauI,
    double const* interIc,
    int const* symptomatic,
    double const* kappa_vals,
    int const* hosp_pop,
    int const* icu_pop,
    double const* lat_locale,
    double const* lon_locale,
    int const* locale_HH,
    int const* HH,
    struct locale const* locale_list,
    int const* severe,
    double const* kappa_cache,

    double* tmp_comm_inf_arr)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < num_infectious) {
        double tmp_comm_inf = 0.0;
        int infec_person; //Counter for infected person.
        double kappa; // #Infectiousness
        double tIc;
        infec_person = infectious[i];
        tIc = Ic;
        if ( intervene[infec_person] > 0 && t>tau[infec_person]+tauI[intervene[infec_person]]) {
            tIc = interIc[intervene[infec_person]];
        }
        kappa = kappa_cache[i];
    
        if (hosp_pop[infec_person]==0) {
            double d; //distance between people.
            // Community transmission //
#if !defined(USE_LOCALE_DISTANCE)
            d = d_distance(lat_locale[j], lon_locale[j], lat_locale[locale_HH[HH[infec_person]]], lon_locale[locale_HH[HH[infec_person]]], 'K');
#else
            d = locale_distance(locale_list[j], locale_list[locale_HH[HH[infec_person]]]);
#endif
            tmp_comm_inf += tIc * calc_community_infect( kappa, omega, severe[infec_person], d, betac_scale);
        }
    
        tmp_comm_inf_arr[i] = tmp_comm_inf;
    }
}

void locale_infectious_step(LoopInvariantData const& lid, int population, int j, int num_households, int num_infectious, double Ic, double t, double dt, double omega, double& out_tmp_comm_inf, int full_kappa, double R0_scale, double betac_scale) {
    if (num_infectious == 0) {
        out_tmp_comm_inf = 0.0;
        return;
    }
    
    thrust::device_vector<double> d_tmp_comm_inf_arr(num_infectious, 0.0);

    // Run kernel
    size_t const THREAD_COUNT = 512;
    size_t const BLOCK_COUNT = (num_infectious + THREAD_COUNT - 1) / THREAD_COUNT;
    locale_infectious_step_kernel<<<BLOCK_COUNT, THREAD_COUNT>>>(
        j,
        num_infectious,
        Ic,
        t,
        dt,
        omega,
        full_kappa,
        R0_scale,
        betac_scale,
        
        lid.infectious.data().get(),
        lid.intervene.data().get(),
        lid.tau.data().get(),
        lid.tauI.data().get(),
        lid.interIc.data().get(),
        lid.symptomatic.data().get(),
        lid.kappa_vals.data().get(),
        lid.hosp_pop.data().get(),
        lid.icu_pop.data().get(),
        lid.lat_locale.data().get(),
        lid.lon_locale.data().get(),
        lid.locale_HH.data().get(),
        lid.HH.data().get(),
        lid.locale_list.data().get(),
        lid.severe.data().get(),
        lid.kappa_cache.data().get(),

        d_tmp_comm_inf_arr.data().get());

    out_tmp_comm_inf = thrust::reduce(d_tmp_comm_inf_arr.begin(), d_tmp_comm_inf_arr.end(), 0.0, thrust::plus<double>());
}

extern "C" void locale_infectious_loop(int num_locale, int population, int num_households, int num_infectious, int* infectious, double Ic, int* intervene, double t, double* tau, double* tauI, double* interIc, int* symptomatic, double dt, double* kappa_vals, int count_kappa_vals, int* hosp_pop, int* icu_pop, double* lat_locale, double* lon_locale, int* locale_HH, int* HH, struct locale* locale_list, double omega, int* severe, int full_kappa, double R0_scale, double betac_scale, double* commun_nom1, double* fd_tot) {

    size_t const num_I = 9;

    // Allocate device and host arrays
    LoopInvariantData lid;
    lid.infectious.assign(infectious, infectious + num_infectious);
    lid.intervene.assign(intervene, intervene + population);
    lid.tau.assign(tau, tau + population);
    lid.tauI.assign(tauI, tauI + num_I);
    lid.interIc.assign(interIc, interIc + num_I);
    lid.symptomatic.assign(symptomatic, symptomatic + population);
    lid.kappa_vals.assign(kappa_vals, kappa_vals + count_kappa_vals);
    lid.hosp_pop.assign(hosp_pop, hosp_pop + population);
    lid.icu_pop.assign(icu_pop, icu_pop + population);
    lid.lat_locale.assign(lat_locale, lat_locale + num_locale);
    lid.lon_locale.assign(lon_locale, lon_locale + num_locale);
    lid.locale_HH.assign(locale_HH, locale_HH + num_households);
    lid.HH.assign(HH, HH + population);
    lid.locale_list.assign(locale_list, locale_list + num_locale);
    lid.severe.assign(severe, severe + population);
    
    thrust::host_vector<double> kappa_cache(num_infectious);
    size_t i;

#ifdef _OPENMP
#pragma omp parallel for private(i) default(shared)
#endif
    for (i = 0; i < num_infectious; ++i) {
        int infec_person = infectious[i];
        kappa_cache[i] = calc_kappa(t, tau[infec_person], symptomatic[infec_person], dt, kappa_vals, hosp_pop[infec_person], icu_pop[infec_person], full_kappa, R0_scale);
    }
    
    lid.kappa_cache = kappa_cache;

	for (int j=0; j<num_locale; j++) {
		double tmp_comm_inf = 0.0;

		locale_infectious_step(lid, population, j, num_households, num_infectious, Ic, t, dt, omega, tmp_comm_inf, full_kappa, R0_scale, betac_scale);

		commun_nom1[j] = tmp_comm_inf / fd_tot[j];
	}
}
