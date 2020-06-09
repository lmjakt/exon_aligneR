#ifndef UTIL_H
#define UTIL_H

#include <unistd.h>  // for size_t

unsigned char *make_complement();
unsigned char *rev_complement(const unsigned char *seq, size_t l, unsigned char *complement);

#endif
