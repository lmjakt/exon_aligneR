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
    return((unsigned char*)0);
  if(l == 0)
    l = strlen((const char*)seq);
  unsigned char *rev_comp = malloc( sizeof(unsigned char) * (l + 1) );
  rev_comp[l] = 0;
  for(int i=0; i < l; ++i){
    rev_comp[l-(i+1)] = complement[ seq[i] ];
  }
  return( rev_comp );
}

// treat '-' as a special character.. 
// note that gap should be negative as we will 'add it on'
void local_score(const char *seq1, const char *seq2, int radius,
		 int gap, int *sub_table, int al_offset, int al_size,
		  double *scores, int al_length)
{
  if(al_length < 1)
    return;

  memset( (void*)scores, 0, sizeof(double) * al_length );

  int score = 0;
  for(int i =0; i < radius && i < al_length; ++i){
    if( seq1[i] == '-' || seq2[i] == '-' )
      score += gap;
    else
      score += sub_table[ (seq1[i] - al_offset) + (seq2[i] - al_offset) * al_size ];
  }    

  for(int i=0; i < al_length; ++i){
    int wbeg = i - radius ;
    int wend = i + radius;

    if(wend < al_length){
      if( seq1[wend] == '-' || seq2[wend] == '-' )
	score += gap;
      else
	score += sub_table[ (seq1[wend] - al_offset) + (seq2[wend] - al_offset) * al_size ];
    }

    if(wbeg >= 0){
      if( seq1[wbeg] == '-' || seq2[wbeg] == '-' )
	score -= gap;
      else
	score -= sub_table[ (seq1[wbeg] - al_offset) + (seq2[wbeg] - al_offset) * al_size ];
    }
    scores[i] = (double)score / (double)( (wend >= al_length ? (al_length-1) : wend) - 
					  (wbeg < 0 ? -1 : wbeg) );
  }
}
