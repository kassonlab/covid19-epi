#ifndef COV_GEO_H_6f1ca88e_624e_4529_8dfa_98eb8cc36d44
#define COV_GEO_H_6f1ca88e_624e_4529_8dfa_98eb8cc36d44

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

struct COV_geo_slice {
	size_t capacity;
	size_t count;
    uint32_t* ids;
    float* lats;
    float* lons;
};

struct COV_geo_band {
	size_t slice_count;
	struct COV_geo_slice* slices;
};

struct COV_geo {
	size_t band_count;
	float lat_min;
	float lat_max;
    float lon_min;
	float lon_max;

	struct COV_geo_band* bands;
};

struct COV_geo_query_context {
	struct COV_geo* geo;
	double lat_cutoff;
	size_t first_band;
	size_t last_band;
	size_t current_band;
	size_t current_slice;
	size_t current_slice_start;
	size_t current_slice_end;
	bool iteration_started;
};

// Input latitude/longitude are in degrees
// Output distance is in kilometers
double COV_geo_distance(double lat1, double lon1, double lat2, double lon2);

struct COV_geo* COV_create_geo(size_t coord_count, uint32_t* ids, float* lat, float* lon,
	size_t band_count, size_t slice_count,
    float lat_min, float lon_min, float lat_max, float lon_max);

void COV_destroy_geo(struct COV_geo* geo);

struct COV_geo_query_context COV_geo_prepare_query(
    struct COV_geo* geo, float lat, float lon, float cutoff_distance);

struct COV_geo_slice* COV_geo_query_next(struct COV_geo_query_context* qc);

#endif /* COV_GEO_H_6f1ca88e_624e_4529_8dfa_98eb8cc36d44 */