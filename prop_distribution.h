#if !defined(__PROP_DISTRIBUTION_H__)
#define __PROP_DISTRIBUTION_H__

long rand_distr_indx(double *dist,
                     long    lb,
                     long    ub);

int generate_inc_distr_vec(double **dist_vec,
                           double  *prob,
                           int      n,
                           char    *name);

#endif /* __PROP_DISTRIBUTION_H__ */
