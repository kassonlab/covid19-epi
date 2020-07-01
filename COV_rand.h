#if !defined __COV_RAND_H__
#define __COV_RAND_H__

/* Wrapper functions for random generator */

extern int use_fixed_seed;

extern long seed;

void   COV_init_rand();

double COV_rand();

#endif /* __COV_RAND_H__ */
