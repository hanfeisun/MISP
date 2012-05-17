#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "kseq.h"
#include "motif.h"
KSEQ_INIT(gzFile, gzread, gzrewind)
#include "debug.h"

#define BITSHIFT 2
#define LOGODD_MIN -256
#define KINDS 4			/* A, T, C, G */
float cg_percent(kseq_t* kseq);
void expected_difference(struct pssm_matrix *pm, float c_p, double* ediff);
void lookahead_filter(int q, kseq_t *kseq, struct pssm_matrix *pm, float c_p, double tol, struct match_doublet* md);


int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	
	FILE* pssm_fp;
	struct pssm_matrix *pssm;
	struct match_doublet *dm;
	
	
	
	pssm_fp = fopen("./database/cistrome.db","r");
	
	pssm = malloc(sizeof(struct pssm_matrix));
	dm = malloc(sizeof(struct match_doublet));
	pssm_reader(pssm_fp, pssm);
	display_pssm(pssm);	
	pssm2logodd(pssm, 0.25);
	display_pssm(pssm);
	
	fclose(pssm_fp);
	/* if (argc == 1) { */
	/* 	fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]); */
	/* 	return 1; */
	/* } */
	/* fp = gzopen(argv[1], "r"); */
	/* seq = kseq_init(fp); */
	
	/* printf("CG percent is %f\n",cg_percent(seq)); */
	/* kseq_rewind(seq); */
	/* lookahead_filter(5, seq, pssm, 0.25, threshold_fromP(pssm, 0.25, 0.1), dm); */
	/* /\* for(;dm->next != NULL; dm=dm->next) *\/ */
	/* /\* 	printf("%d,%f\n",dm->position,dm->score); *\/ */

	/* kseq_destroy(seq); */
	/* gzclose(fp); */
	return 0;
}

float cg_percent(kseq_t* kseq)
{
	int l;
	unsigned long a, c, other;
	a = 0;
	c = 0;
	other = 0;
	while ((l = kseq_read(kseq)) >= 0) 
		bg_counter(&a, &c, &other, kseq->seq.s);

	printf("A or T: %ld\n", a);
	printf("C or G: %ld\n", c);
	printf("Others: %ld\n", other);
	
	if (a+c<50) {
		printf("Your input sequence has too few lines, the bg score will be set to 0.5\n");
		return 0.5;
	} else  return ((double)a)/(a+c);
}


void lookahead_filter(int q, kseq_t *kseq, struct pssm_matrix *pm, float c_p, double tol, struct match_doublet* md)
{
	int i, j;
	int window_pos;
	int *order;
	int *order_y;
	double *good;
	double *scores;
	double tmp, tmp_max;
	int bufsize;
	short *seq_code;
	uint_fast32_t *sA;
	uint_fast32_t code;
	const uint_fast32_t size = 1 << (BITSHIFT * q);;

	
	assert(pm->kinds == KINDS);	/* Only for DNA */
	/* Find the window */
	if (q > pm->len) {
		lookahead_filter(pm->len, kseq, pm, c_p, tol, md);
		return;
	}
	sA = calloc(q, sizeof(uint_fast32_t));	
	good = malloc(sizeof(double) * pm->len);
	tmp = 0;
	expected_difference(pm, c_p, good);
	for (i = 0; i < pm->len; ++i)
		tmp += good[i];
	tmp_max = tmp;
	window_pos = 0;

	for (i = 0; i < pm->len - q; ++i) {
		tmp -= good[i];
		tmp += good[i+q];
		if (tmp > tmp_max) {
			tmp_max = tmp;
			window_pos = i + 1;
		}
	}

	/* Arrange matrix indeces remained */
	order = calloc(pm->len - q, sizeof(int));
	for (i = 0; i < window_pos; ++i)
		order[i] = i;
	for (i = window_pos + q; i < pm->len; ++i)
		order[i-q] = i;
	qsort(order, pm->len - q, sizeof(int), compare);

	/* Lookahead array for indeces outside the window */
	free(good);
	good = malloc(sizeof(double) * (pm->len - q));
	for (j = pm->len - q; j > 0; --j) {
		tmp_max = LOGODD_MIN;
		for (i = 0; i < pm->kinds; ++i) {
			if (tmp_max < pm->score[i][order[j]])
				tmp_max = pm->score[i][order[j]];
		}
		good[j-1] = good[j] + tmp_max;
	}

	/* Scores in the window for all possible strings */
	
	scores = malloc(sizeof(double) * (size));
	tmp = 0;
	for (j = 0; j < q; ++j)
		tmp = tmp + pm->score[0][window_pos+j];
	scores[0] = tmp;
	/* printf("%f\n",scores[0]); */
	j = 0;
	code = 0;
	
	while(1) {
		++code;
		if (code >= size)
			break;

		scores[code] = scores[code - 1];
		j = q - 1;	/* j is pos now */
		while (j >= 0) {
			if (sA[j] <  KINDS - 1) {
				scores[code] +=				\
					pm->score[sA[j] + 1][j + window_pos] \
					-				\
					pm->score[sA[j]][j + window_pos];
				++sA[j];
				break;
			} else {
				scores[code] +=				\
					pm->score[0][j + window_pos]	\
					-				\
					pm->score[KINDS - 1][j + window_pos];
				sA[j] = 0;
				--j;
			}

		}
		/* printf("%d's score is %f\n", code, scores[code]); */
	}
	
	/* Actual scanning */
	tmp = 0;
	tmp_max = tol - good[0];


	while ((bufsize = kseq_read(kseq)) >= 0) {
		seq_code = malloc(sizeof(short) * (bufsize+1));

		/* printf("base2code\n"); */
		base2code(kseq->seq.s, seq_code);

		/* for (i = 0; seq_code[i] != -1; ++i) */
		/* 	printf("%d",seq_code[i]); */
		/* printf("\n%s\n",kseq->seq.s); */
	
		code = 0;


		for (i = window_pos; i < window_pos + q - 1; ++i) {
			code = (code << BITSHIFT) + seq_code[i];
		}
	

		for(i = 0; seq_code[pm->len - 1]!=-1; ++seq_code) {
			++i;
			code = ((code << BITSHIFT) + seq_code[window_pos + q -1]) & (size - 1);
			/* printf("now code is %d\n",code); */

			if (scores[code] >= tmp_max) {
				tmp = scores[code];
				order_y = order;
				for (j = 0; j < pm->len - q; ++j) {
					if (tmp + good[j] < tol)
						break;
					tmp += pm->score[seq_code[*order_y]][*order_y];
					++order_y;
				}
				if (tmp >= tol) {
					md->position = i;
					md->score = tmp;
					md->next = malloc(sizeof(struct match_doublet));
					md = md->next;
				}
			
			}
		
		}
	}
	
	md->next = NULL;
}

int compare (const void * a, const void * b)
{
	return (*(int*)a - *(int*)b);
}

void expected_difference(struct pssm_matrix *pm, float c_p, double* ediff)
/* Please allocate enough memory for ediff before call me! */
{
	int i,j;
	float max;

	for (i = 0; i < pm->len; ++i) {
		max = LOGODD_MIN;
		for (j = 0; j < pm->kinds; ++j) {
			if (max < pm->score[j][i])
				max = pm->score[j][i];
		}
 
		ediff[i] = max;

		for (j = 0; j < pm->kinds; ++j) {
			if (j == PWM_BASE_C || j == PWM_BASE_G)
				ediff[i] -= c_p * pm->score[j][i];
			else if (j == PWM_BASE_A || j == PWM_BASE_T)
				ediff[i] -= (0.5 - c_p) * pm->score[j][i];
			else {
				perror("Error");
				exit(1);
			}
		}
	}
}

	
	
