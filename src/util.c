#include "util.h"
#include <ctype.h>  // for toupper
#include <stdlib.h> // for malloc
#include <string.h> // for strlen

unsigned char *make_complement(){
  unsigned char *complement = malloc(sizeof(unsigned char) * 256);
  for(size_t c=0; c < 256; c++)
    complement[c] = c;
  // from: https://droog.gs.washington.edu/parc/images/iupac.html
  // IUPAC ambiguity codes and complements
  unsigned char nucs[] = {'a', 'c', 'g', 't', 'm', 'r', 'w', 's', 'y', 'k', 'v', 'h', 'd', 'b'}; 
  unsigned char comp[] = {'t', 'g', 'c', 'a', 'k', 'y', 'w', 's', 'r', 'm', 'b', 'd', 'h', 'v'};
  int l = 14;
  for(int i=0; i < l; ++i){
    complement[ nucs[i] ] = comp[i];
    complement[ toupper(nucs[i]) ] = toupper( comp[i] );
  }
  return(complement);
}

unsigned char *rev_complement(const unsigned char *seq, size_t l, unsigned char *complement){
  if(seq == 0 || complement == 0)
    return((const char*)0);
  if(l == 0)
    l = strlen((const char*)seq);
  unsigned char *rev_comp = malloc( sizeof(unsigned char) * (l + 1) );
  rev_comp[l] = 0;
  for(int i=0; i < l; ++i){
    rev_comp[l-(i+1)] = complement[ seq[i] ];
  }
  return( rev_comp );
}
