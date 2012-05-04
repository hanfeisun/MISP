#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "motif.h"
KSEQ_INIT(gzFile, gzread)
#include "debug.h"

float cg_percent(kseq_t* kseq);
float* lookahead_filter(int q, kseq_t *kseq, struct pwm_matrix **mpm, float cg_p, double tol);



int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	
	FILE* pwm_fp;
	struct pwm_matrix *pwm;
	struct pssm_matrix *pssm;
	
	
	pwm_fp = fopen("MA0064.pfm","r");
	
	pwm = malloc(sizeof(struct pwm_matrix));
	pssm = malloc(sizeof(struct pssm_matrix));
	pwm_reader(pwm_fp,pwm);
	display_pwm(pwm);
	counts2Logfodds(pwm, pssm, 0.25, PSUEDO_E);
	display_pssm(pssm);
	printf("threshold is %f\n",threshold_fromP(pssm, 0.25, 0.99));
	printf("threshold is %f\n",threshold_fromP(pssm, 0.25, 0.1));
	printf("threshold is %f\n",threshold_fromP(pssm, 0.25, 0.00001));	
	
	fclose(pwm_fp);
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	printf("CG percent is %f\n",cg_percent(seq));
	kseq_rewind(seq);
	printf("CG percent is %f\n",cg_percent(seq));
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

	printf("A or T: %ld\n", a);
	printf("C or G: %ld\n", c);
	printf("Others: %ld\n", other);
	
	if (a+c<50) {
		printf("Your input sequence has too few lines, the bg score will be set to 0.5\n");
		return 0.5;
	} else  return ((double)a)/(a+c);
}


