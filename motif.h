#ifndef MIS_MOTIF_H
#define MIS_MOTIF_H
#define PWM_MAX_COL 50
#define PVAL_DP_MULTIPLIER 1000.0
#define PWM_BASE_A 0
#define PWM_BASE_C 1
#define PWM_BASE_G 2
#define PWM_BASE_T 3
#define PSUEDO_E 0.1
#include "kseq.h"
#include <zlib.h>
#include <stdio.h>



struct pwm_matrix 
{
	int **weight;
	/* weight[kinds][position] */
	int len;
	short kinds;
	/* char *comment; */
};

struct pssm_matrix 
{
	float **score;
	/* weight[kinds][position] */
	int len;
	short kinds;
	/* char *comment; */
};

struct match_pos 
{
	int position;
	float score;
};

int bg_counter(unsigned long *acnt, unsigned long *ccnt, unsigned long *other, char *seq);

int pwm_reader(FILE *fp, struct pwm_matrix *pm);

void display_pwm(struct pwm_matrix *pm);

void display_pssm(struct pssm_matrix *pm);

// Calculates a threshold for a scoring matrix from a given p value
double threshold_fromP (struct pssm_matrix *pm, double c_p, double p);

inline int init_array(double *array, size_t size, double value);

int counts2Logfodds(struct pwm_matrix *pm_in, struct pssm_matrix *pm_out,  double c_p, double ps);

#endif
