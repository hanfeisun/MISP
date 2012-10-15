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
void lookahead_filter(int q, kseq_t *kseq, struct pssm_matrix *pm, float c_p, double tol, struct match_doublet *md, FILE *output, int mini);


uint_fast32_t flip_reverse2(uint_fast32_t value)
{
	int i;
	uint32_t new_value = 0;
	for (i = 0; i < 16; i++)
	{
		new_value <<= 2;            
		new_value |= (value & 0x3);
		value >>= 2;
	}
	return ((~new_value)>>22) & 0x3ff; /* only ok when q is 5 (22 == 32 - 2*5)*/
}



int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	FILE* pssm_fp;
	FILE* output;
	struct pssm_matrix *pssm;
	char of_name[2000];
	struct match_doublet *dm;
	float p_cg;
	float p_value;

	pssm = malloc(sizeof(struct pssm_matrix));
	dm = malloc(sizeof(struct match_doublet));

	if (argc < 5) {
		fprintf(stderr, "Usage: %s <in.seq> <in.db> <p-value> <motif-id> <output-prefix>\n", argv[0]);
		return 1;
	}
	
	pssm_fp = fopen(argv[2],"r");
	pssm_reader(pssm_fp, pssm);
	fclose(pssm_fp);
	
	fp = gzopen(argv[1], "r");
	p_value = atof(argv[3]);

	seq = kseq_init(fp);
	p_cg=cg_percent(seq);
	pssm2logodd(pssm, p_cg);
	
	if (strcasecmp("all", argv[4]) == 0) {
		sprintf(of_name,"%s_all",argv[5]);
		output = fopen(of_name,"w");
		fprintf(output, "# Parameter List:\n");
		fprintf(output, "# Input sequence: %s\n",argv[1]);
		fprintf(output, "# Output path: %s\n",argv[5]);
		fprintf(output, "# P value: %.3f\n",p_value);	
		for(; pssm->next != NULL; pssm = pssm->next) {
			kseq_rewind(seq);
			lookahead_filter(5, seq, pssm, p_cg/2.0, threshold_fromP(pssm, p_cg/2.0, p_value), dm, output,1);
		}
		fclose(output);
		
	} else {
		for(; pssm->next != NULL; pssm = pssm->next) {
			if (strncasecmp(pssm->name, argv[4], strlen(pssm->name) - 1) == 0)
				break;
		}

		if (strncasecmp(pssm->name, argv[4], strlen(pssm->name) -1) != 0) {
			printf("Can't find motif id %s", argv[3]);
			return 0;
		}
		kseq_rewind(seq);
		pssm->name[strlen(pssm->name) - 1] = '\0';
		sprintf(of_name,"%s_%s",argv[5],pssm->name);
		output = fopen(of_name,"w");
		fprintf(output, "# Parameter List:\n");
		fprintf(output, "# Input sequence: %s\n",argv[1]);
		fprintf(output, "# Output path: %s\n",argv[5]);
		fprintf(output, "# P value: %.3f\n",p_value);	
		lookahead_filter(5, seq, pssm, p_cg/2.0, threshold_fromP(pssm, p_cg/2.0, p_value), dm, output, 0);
		fclose(output);
	}
	kseq_destroy(seq);
	gzclose(fp);
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
	printf("Calculating Background..\n");
	printf("A or T: %ld\n", a);
	printf("C or G: %ld\n", c);
	printf("Others: %ld\n", other);
	if (a+c<50) {
		printf("Your input sequence has too few lines, the bg score will be set to 0.5\n");
		return 0.5;
	} else  return ((double)c)/(a+c);
}


void lookahead_filter(int q, kseq_t *kseq, struct pssm_matrix *pm, float c_p, double tol, struct match_doublet* md, FILE *output, int mini)
{
	int i, j;
	int window_pos;
	struct order_s **order;
	struct order_s **order_y;
	double *good;
	double *scores;
	int pos_max;
	double hit_max;
	double tmp, tmp_max;
	int bufsize;
	int has_hit;
	short *seq_code;
	short *seq_code_gc;
	uint_fast32_t *sA;
	uint_fast32_t code;
	const uint_fast32_t size = 1 << (BITSHIFT * q);;
	int compare (const void * a, const void * b);
	
	assert(pm->kinds == KINDS);	/* Only for DNA */
	/* Find the window */
	if (q > pm->len) {
		return lookahead_filter(pm->len, kseq, pm, c_p, tol, md, output,mini);
	}
	sA = calloc(q, sizeof(uint_fast32_t));	
	good = malloc(sizeof(double) * pm->len);
	tmp = 0;
	expected_difference(pm, c_p, good);
	for (i = 0; i < q; ++i)
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
	/* NOTICE: start from 1 !!*/
	order = (struct order_s **) calloc(pm->len - q, sizeof(struct order_s*));
	
	for (i = 0; i < window_pos; ++i) {
		order[i] = (struct order_s *) malloc(sizeof(struct order_s));
		order[i]->pos = i+1;
		order[i]->good = good[i];
	}
	for (i = window_pos + q; i < pm->len; ++i)
	{
		order[i-q] = (struct order_s *) malloc(sizeof(struct order_s));
		order[i-q]->pos = i+1;
		order[i-q]->good = good[i];
	}
	qsort(order, pm->len - q, sizeof(struct order_s *), compare);
	/* Lookahead array for indeces outside the window */
	good = malloc(sizeof(double) * (pm->len - q));
	for (j = pm->len - q; j > 0; --j) {
		tmp_max = LOGODD_MIN;
		for (i = 0; i < pm->kinds; ++i) {
			if (tmp_max < pm->score[i][order[j-1]->pos])
				tmp_max = pm->score[i][order[j-1]->pos];
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
	

	if (mini == 0) {
		fprintf(output, "# CG percent: %.2f\n", c_p*2);
		fprintf(output, "# tolerance: %.2f\n", tol);
		fprintf(output, "# factor ID: %s\n", pm->name);
		fprintf(output, "sequence_name\tsequence_length\thits_score\thits_position\tsequence\n");
	} else {
		fprintf(output, "\n\n# factor:%s# tolerance %.2f\n", pm->name, tol);
	}
	
	while ((bufsize = kseq_read(kseq)) >= 0) {
		if (mini == 0) 
			fprintf(output,"%s\t%d\t", kseq->name.s, bufsize);
			
		if (bufsize < pm->len) {
			if (mini == 0)
				fprintf(output, "*\t*\t*\n");
			continue;
		}

		seq_code = malloc(sizeof(short) * (bufsize+1));
		seq_code_gc = seq_code;

		base2code(kseq->seq.s, seq_code);
	
		code = 0;

		int N_max_in_window=0;
		for (i = window_pos; i < window_pos + q - 1; ++i) {
			if (seq_code[i] == 4)
				N_max_in_window = i;
			code = (code << BITSHIFT) + seq_code[i];
		}

		hit_max = -10;
		has_hit = 0;
		pos_max = -1;
		int positive_end = 0;
		int negative_start = 0;
		for(i = 0; seq_code[window_pos + q - 1]!=-1; ++i && ++seq_code) {
			/* i starts from 0, but pos_max also starts from 0 */
			code = ((code << BITSHIFT) + seq_code[window_pos + q -1]) & (size - 1);
			if (N_max_in_window != 0) {
				N_max_in_window -= 1; /* skip  */
				break;
			}
			if (seq_code[window_pos + q -1] == 4) {
				N_max_in_window = q - 1;
				break;
			}
			
			  
			/* printf("now code is %d\n",code); */
			if (seq_code[(pm->len)-1] == -1) 
				positive_end = 1;
			
			
			if (!positive_end && ((tmp = scores[code]) >= tmp_max))  { /* tmp max is the lower bound */
				order_y = order;
				for (j = 0; j < pm->len - q; ++j) {
					if (tmp + good[j] < tol)
						break;
					if (seq_code[order_y[0]->pos - 1] == 4) {
						tmp = tol - 1; /* abandon this sequence if N exists */
					  /* For masked */
						break;
					}
					
					tmp += pm->score[seq_code[order_y[0]->pos - 1]][order_y[0]->pos-1];
					++order_y;

				}
				if (tmp >= tol) {
					tmp -= tol;
					if (tmp > hit_max) {
						hit_max = tmp;
						pos_max = i;
					}
					has_hit = 1;

				}

			}
			if (i >= pm->len - q)
				negative_start = 1;
			/* pm->len >= window_pos + q */
			if (negative_start && ((tmp = scores[flip_reverse2(code)]) >= tmp_max)) { /* reverse strand */
				order_y = order;
				for (j = 0; j < pm->len - q; ++j) {
					if (tmp + good[j] < tol)
						break;
					if (seq_code[q - order_y[0]->pos ] == 4) { /* for masked Basepair */
						tmp = tol - 1;
						break;
					}
					tmp += pm->score[3 - seq_code[q + window_pos - order_y[0]->pos]][order_y[0]->pos - 1];
					++order_y;
				}
				if (tmp >= tol) {
					tmp -= tol;
					if (tmp > hit_max) {
						hit_max = tmp;
						pos_max = -i;
					}
					has_hit = 1;
				}
			}
				
				
		}

		free(seq_code_gc);
		if (mini == 0) {
		if (has_hit){
				char temp_char;
				int e_idx, s_idx; /* index of end and start position */
				if (pos_max>0)
					e_idx = pos_max + pm->len;
				else
					e_idx = - pos_max + q + window_pos;
				s_idx = e_idx - pm->len;
				temp_char = *(kseq->seq.s + e_idx);
				*(kseq->seq.s + e_idx) = '\0';
				if (pos_max>0) 
					fprintf(output, "\t%.2f\t%d\t%s\n",hit_max, pos_max, kseq->seq.s + s_idx);
				else 
					fprintf(output, "\t%.2f\t-%d\t%s\n",hit_max, e_idx , kseq->seq.s + s_idx);
				*(kseq->seq.s + e_idx) = temp_char;
			}
			else
				fprintf(output, "*\t0\t*\n");
		} else {
			if (has_hit)
				fprintf(output, "%.2f(%d),", hit_max, pos_max);
			else
				fprintf(output, "0,");
		}
		
	}
}



int compare (const void * a, const void * b)
{
	if (((*(const struct order_s **)a)->good - (*(const struct order_s **)b)->good) > 0)
		return 1;
	else
		return -1;
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

	
	
