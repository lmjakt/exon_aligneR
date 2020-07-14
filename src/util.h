#ifndef UTIL_H
#define UTIL_H

#include <unistd.h>  // for size_t

unsigned char *make_complement();
unsigned char *rev_complement(const unsigned char *seq, size_t l, unsigned char *complement);

void local_score(const char *seq1, const char *seq2, int radius,
		 int gap, int *sub_table, int al_offset, int al_size,
		 double *scores, int al_length);

#endif
