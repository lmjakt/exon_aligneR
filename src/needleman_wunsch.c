#include <string.h>
#include <stdlib.h>
#include <string.h>
#include "needleman_wunsch.h"
#include <Rinternals.h>

// A very basic function that takes
// two sequences
// gap insertion penalty
// gap extension penalty
// a substitution table
// an offset (the smallest value of the alphabet)
// the alphabet offset
// the alphabet size
// a logical value indicating whether terminal gaps in the
// shorter sequence should be penalised.
//
// The calling function should also assign memory for the
// score and pointer tables. The function does not return anything
// but instead fills in these tables which can then be used by
// some other functions for extracting the desired alignment.
// 
// the substitution matrix should be of size:
//     size * size
// and row and column indices should be calculated by:
//
//  r - offset
//  where r is the char value of the sequence.
// This function does not check that this will give a valid
// column number.
//
// the matrix is row-major as is standard for R.
//
// The basic idea here is to insert special characters at intron positions
// And to give these very high alignment scores.
// This should favour alignments that are exon-aligned.
//
// Conceptually this could be used for other purposes as well.


int which_max_i(int *v, int l){
  int max_i = 0;
  double max = v[0];
  for(int i=1; i < l; i++){
    if(v[i] > max){
      max = v[i];
      max_i = i;
    }
  }
  return(max_i);
}

int which_max_d(double *v, int l){
  int max_i = 0;
  double max = v[0];
  for(int i=1; i < l; i++){
    if(v[i] > max){
      max = v[i];
      max_i = i;
    }
  }
  return(max_i);
}


// To enforce the matrices; always use this
int m_offset(int row, int column, int height){
  return( column * height + row );
}

void init_scores( int *scores, int height, int width, char gap_i, char gap_e, char tgaps_free ){
  //  Rprintf("init scores dims: %d, $d\n", height, width );
  memset( (void*)scores, 0, sizeof(int) * height * width );
  if( width >= height || !(tgaps_free) ){
    scores[1] = gap_i;
    for(int i=2; i < height; ++i)
      scores[i] = scores[i-1] + gap_e;
  }
  if( height >= width || !(tgaps_free)){
    scores[ height ] = gap_i;
    for(int i=2; i < width; ++i)
      scores[ i * height ] = scores[ (i-1) * height ] + gap_e;
  }
}

void init_ptr( int *ptr, int height, int width, const int left, const int up ){
  memset( (void*)ptr, 0, sizeof(int) * height * width );
  //  Rprintf("init pointers dims: %d, $d\n", height, width );
  for(int i=1; i < height; ++i)
    ptr[i] = up;
  for(int i=1; i < width; ++i)
    ptr[ i * height ] = left;
}

// a is represented by rows, b, by columns, to be consistent with R
void needleman_wunsch( const unsigned char *a, const unsigned char *b, int a_l, int b_l,
		       int gap_i, int gap_e,
		       int *sub_table, int al_offset, int al_size,
		       int tgaps_free,
		       int *score_table, int *ptr_table)
{
  int height = a_l + 1;
  int width = b_l + 1;
  const int left = 1;
  const int up = 2;
  /* Rprintf("length of a: %d or %d\n", a_l, strlen((const char*)a)); */
  /* Rprintf("length of b: %d or %d\n", b_l, strlen((const char*)b)); */
  //  Rprintf("Calling init_scores %p  %d, %d   %d, %d\n", score_table, height, width, gap_i, gap_e);
  init_scores( score_table, height, width, gap_i, gap_e, tgaps_free );
  //  Rprintf("Calling int_ptr\n");
  init_ptr( ptr_table, height, width, left, up );
  
  // whether we allow terminal gaps at the beginning or not.
  int a_tgap_free = tgaps_free && a_l > b_l;
  int b_tgap_free = tgaps_free && b_l > a_l;

  //  Rprintf("a_tgap_free: %d b_tgap_free: %d\n", a_tgap_free, b_tgap_free);
  
  for(int row=1; row < height; ++row){
    //    Rprintf("\nrow: %d\n", row);
    for(int column=1; column < width; ++column){
      int o = column * height + row;
      int o_up = o - 1;
      int o_left = o - height;
      int o_diag = o - (height + 1);
      int m_score_o = (a[row-1] - al_offset) + (b[column-1] - al_offset) * al_size;
      //      Rprintf("\t%c:%c %d, %d", a[row-1], b[column-1],  m_score_o, sub_table[ m_score_o ]);
      char m_score = sub_table[ m_score_o ];
      
      char left_penalty = (b_tgap_free && row == a_l) ? 0 : (ptr_table[o_left] == left ? gap_e : gap_i);
      char up_penalty = (a_tgap_free && column == b_l) ? 0 : (ptr_table[o_up] == up ? gap_e : gap_i);
      
      int scores[3];
      // scores are, left, up, diagonal
      scores[0] = score_table[o_left] + left_penalty;
      scores[1] = score_table[o_up] + up_penalty;
      scores[2] = score_table[o_diag] + m_score;
      int max_i = which_max_i( scores, 3 );
      score_table[o] = scores[ max_i ];
      ptr_table[o] = max_i + 1;
    }
  }
  //  Rprintf("\nGot to the bottom of that\n");
  // and that should basically be it.. apart from the initiation functions.
  
}

void extract_nm_alignment(int *pointers, int height, int width, const char *a, const char *b,
			  char **a_a, char **b_a){
  int al_length = 0;
  int row = height-1; // max.row;
  int column = width-1; // max.column;
  while(row > 0 || column > 0){
    int o = m_offset( row, column, height );
    if(pointers[o] == 0)
      break;
    row = (pointers[o] &  2) ? row - 1 : row;
    column = (pointers[o] & 1) ? column - 1 : column;
    al_length++;
  }
  *a_a = malloc(sizeof(char) * (al_length + 1));
  *b_a = malloc(sizeof(char) * (al_length + 1));
  (*a_a)[al_length] = 0;
  (*b_a)[al_length] = 0;
  row = height-1; //max.row;
  column = width-1; // max.column;
  while(row > 0 || column > 0){
    al_length--;
    int o = m_offset( row, column, height );
    if(pointers[o] == 0 || al_length < 0)
      break;
    // This asks the same question twice and is hence bad. But the
    // only alternative I can think of makes use of two if-else
    // constructs and that too is ugly.
    (*a_a)[al_length] = (pointers[o] & 2) ? a[row - 1] : '-';
    (*b_a)[al_length] = (pointers[o] & 1) ? b[column - 1] : '-';
    row = (pointers[o] &  2) ? row - 1 : row;
    column = (pointers[o] & 1) ? column - 1 : column;
  }
}

// finds the positions of occurences of the character c in word
void char_at(const char *word, char c, int **pos, int *pos_l){
  char *w;
  w = (char*)word;
  *pos = 0;
  *pos_l = 0;
  while( (*w) ){
    if( (*w) == c )
      (*pos_l)++;
    ++w;
  }
  if(*pos_l){
    int pos_i = 0;
    int i=0;
    (*pos) = malloc( sizeof(int) * (*pos_l));
    w = (char*)word;
    while( (*w) ){
      if((*w) == c){
	(*pos)[pos_i] = i;
	++pos_i;
      }
      ++i;
      ++w;
    }
  }
}

// pos1 and pos2 are sorted vectors giving the positions of some
// residue in aligned sequences. l1 and l2 are the respective lengths.
// This function returns a table of indices which are aligned to each other
int *aligned_i(int *pos1, int* pos2, int l1, int l2, int *nrow){
  int i1=0;
  int i2=0;
  *nrow = 0;
  int *i_table = 0;
  while(i1 < l1 && i2 < l2){
    if( pos1[i1] == pos2[i2] )
      (*nrow)++;
    if( pos1[i1] <= pos2[i2] )
      i1++;
    else
      i2++;
  }
  if(*nrow){
    i_table = malloc(sizeof(int) * 2 * (*nrow));
    i1 = i2 = 0;
    int row=0;
    while(i1 < l1 && i2 < l2){
      if( pos1[i1] == pos2[i2] ){
	i_table[row] = i1;
	i_table[(*nrow) + row] = i2;
	++row;
      }
      if( pos1[i1] <= pos2[i2] )
	i1++;
      else
	i2++;
    }
  }
  return(i_table);
}
 
