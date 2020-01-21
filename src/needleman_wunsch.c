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

/* int which_max_i(int *v, int l){ */
/*   int max_i = 0; */
/*   double max = v[0]; */
/*   for(int i=1; i < l; i++){ */
/*     if(v[i] > max){ */
/*       max = v[i]; */
/*       max_i = i; */
/*     } */
/*   } */
/*   return(max_i); */
/* } */

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

// Does not init whole table
void init_scores( int *scores, int height, int width, char gap_i, char gap_e, char tgaps_free ){
  //  Rprintf("init scores dims: %d, $d\n", height, width );
  memset( (void*)scores, 0, sizeof(int) * height * width );
  scores[0] = 0;
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

// Does not init the whole table
void init_ptr( int *ptr, int height, int width, const int left, const int up ){
  memset( (void*)ptr, 0, sizeof(int) * height * width );
  ptr[0] = 0;
  //  Rprintf("init pointers dims: %d, $d\n", height, width );
  for(int i=1; i < height; ++i)
    ptr[i] = up;
  for(int i=1; i < width; ++i)
    ptr[ i * height ] = left;
}

struct align_stats align_stats_init(){
  struct align_stats as;
  as.al_length = 0;
  as.al_n = 0;
  as.a_gap_i = 0;
  as.b_gap_i = 0;
  as.a_gap = 0;
  as.b_gap = 0;
  as.a_gap_l = 0;
  as.a_gap_r = 0;
  as.b_gap_l = 0;
  as.b_gap_r = 0;
  as.match_n = 0;
  as.mismatch_n = 0;
  as.transition = 0;
  as.transversion = 0;
  as.A = 0;
  as.C = 0;
  as.G = 0;
  as.T = 0;
  // better or worse as:
  //  memset( (void*)&as, 0, sizeof(as) );
  return(as);
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

  init_scores( score_table, height, width, gap_i, gap_e, tgaps_free );
  init_ptr( ptr_table, height, width, left, up );
  
  // whether we allow terminal gaps at the beginning or not.
  int a_tgap_free = tgaps_free && a_l > b_l;
  int b_tgap_free = tgaps_free && b_l > a_l;

  for(int row=1; row < height; ++row){
    for(int column=1; column < width; ++column){
      int o = column * height + row;
      int o_up = o - 1;
      int o_left = o - height;
      int o_diag = o - (height + 1);
      int m_score_o = (a[row-1] - al_offset) + (b[column-1] - al_offset) * al_size;
      int m_score = sub_table[ m_score_o ];
      
      int left_penalty = (b_tgap_free && row == a_l) ? 0 : (ptr_table[o_left] == left ? gap_e : gap_i);
      int up_penalty = (a_tgap_free && column == b_l) ? 0 : (ptr_table[o_up] == up ? gap_e : gap_i);
      
      // This is marginally faster than defining an array of value. Up to 20%.
      score_table[o] = score_table[o_left] + left_penalty;
      ptr_table[o] = 1;
      if( score_table[o_up] + up_penalty > score_table[o] ){
	score_table[o] = score_table[o_up] + up_penalty;
	ptr_table[o] = 2;
      }
      if( score_table[o_diag] + m_score > score_table[o] ){
	score_table[o] = score_table[o_diag] + m_score;
	ptr_table[o] = 3;
      }
      
    }
  }
}

struct align_stats extract_nm_alignment(int *pointers, int height, int width, const unsigned char *a, const unsigned char *b,
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
  // Let us count gaps, gap insertions, matches, mismatches
  char last_a=0, last_b=0;  //
  char ac=0, bc=0;
  struct align_stats stats = align_stats_init();
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
    // The alignment stats
    ac = (*a_a)[al_length];
    bc = (*b_a)[al_length];
    stats.al_length++;
    stats.al_n += ( ac != '-' && bc != '-' ) ? 1 : 0;
    if( ac == '-' ){
      stats.a_gap++;
      stats.a_gap_i += ( last_a == '-' ? 0 : 1 );
      stats.a_gap_l += (row == 0) ? 1 : 0;
      stats.a_gap_r += (row == height-1) ? 1 : 0;
    }
    if( bc == '-' ){
      stats.b_gap++;
      stats.b_gap_i += ( last_b == '-' ? 0 : 1 );
      stats.b_gap_l += (column == 0) ? 1 : 0;
      stats.b_gap_r += (column == width - 1) ? 1 : 0;
    }
    int m_type = mut_type(ac, bc, &stats.A, &stats.C, &stats.G, &stats.T);
    if( ac == bc ){
      stats.match_n++;
    }else{
      stats.mismatch_n += (ac == '-' || bc == '-') ? 0 : 1;
      stats.transition += (m_type == -1) ? 1 : 0;
      stats.transversion += (m_type == 1) ? 1 : 0;
    }
    // and.. 
    last_a = ac;
    last_b = bc;
  }
  return(stats);
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

enum nuc_combination {
  AA = ('A' << 8) | 'A',
  AC = ('A' << 8) | 'C',
  AG = ('A' << 8) | 'G',
  AT = ('A' << 8) | 'T',

  CA = ('C' << 8) | 'A',
  CC = ('C' << 8) | 'C',
  CG = ('C' << 8) | 'G',
  CT = ('C' << 8) | 'T',

  GA = ('G' << 8) | 'A',
  GC = ('G' << 8) | 'C',
  GG = ('G' << 8) | 'G',
  GT = ('G' << 8) | 'T',

  TA = ('T' << 8) | 'A',
  TC = ('T' << 8) | 'C',
  TG = ('T' << 8) | 'G',
  TT = ('T' << 8) | 'T',
};


// this is a bit ugly, but should be fast as it will make use of
// single flat switch statement.
int mut_type(char a, char b, int *A, int *C, int *G, int *T){
  int mask = ~0x20;  // implements toupper.
  int n_combo = ((a & mask) << 8) | (b & mask);
  int m_type = 0;
  switch(n_combo){
  case AA:
    (*A) += 2;
    break;
  case AC:
    (*A)++;
    (*C)++;
    m_type = 1;
    break;
  case AG:
    m_type = -1;
    (*A)++;
    (*G)++;
    break;
  case AT:
    (*A)++;
    (*T)++;
    m_type = 1;
    break;

  case CA:
    (*C)++;
    (*A)++;
    m_type = 1;
    break;
  case CC:
    (*C) += 2;
    break;
  case CG:
    (*C)++;
    (*G)++;
    m_type = 1;
    break;
  case CT:
    (*C)++;
    (*T)++;
    m_type = -1;
    break;

  case GA:
    (*G)++;
    (*A)++;
    m_type = -1;
    break;
  case GC:
    (*G)++;
    (*C)++;
    m_type = 1;
    break;
  case GG:
    (*G) += 2;
    break;
  case GT:
    (*G)++;
    (*T)++;
    m_type = 1;
    break;

  case TA:
    (*T)++;
    (*A)++;
    m_type = 1;
    break;
  case TC:
    (*T)++;
    (*C)++;
    m_type = -1;
    break;
  case TG:
    (*T)++;
    (*G)++;
    m_type = 1;
    break;
  case TT:
    (*T) += 2;
    break;
  default:
    m_type = 0;
  }
  return(m_type);
}

// these functions have too many arguments. It would be better to define suitable structs
// and quite possibly return them.
void smith_waterman( const unsigned char *a, const unsigned char *b, int a_l, int b_l,
		     int gap_i, int gap_e,
		     int *sub_table, int al_offset, int al_size,
		     int *score_table, int *ptr_table,
		     int *max_score, int *max_row, int *max_column)
{
  int height = a_l + 1;
  int width = b_l + 1;
  const int left = 1;
  const int up = 2;
  const int diag = 3;

  // for a smith waterman the initial score table is all 0s as we can start
  // an alignment from anywhere. 
  memset( (void*)score_table, 0, sizeof(int) * height * width );
  memset( (void*)ptr_table, 0, sizeof(int) * height * width );

  *max_score = 0;
  *max_row = 0;
  *max_column = 0;
  
  for(int row=1; row < height; ++row){
    for(int column=1; column < width; ++column){
      int o = column * height + row;
      int o_up = o - 1;
      int o_left = o - height;
      int o_diag = o - (height + 1);
      int m_score_o = (a[row-1] - al_offset) + (b[column-1] - al_offset) * al_size;

      int left_score = score_table[o_left] + (ptr_table[o_left] == left) ? gap_e : gap_i;
      int up_score = score_table[o_up] + (ptr_table[o_up] == up);
      int diag_score = score_table[o_diag] + sub_table[ m_score_o ];

      // we could probably do this by directly calculating the alternative scores
      // rather than
      if(left_score > score_table[o]){
	score_table[o] = left_score;
	ptr_table[o] = left;
      }
      if(up_score > score_table[o]){
	score_table[o] = up_score;
	ptr_table[o] = up;
      }
      if(diag_score > score_table[o]){
	score_table[o] = diag_score;
	ptr_table[o] = diag;
      }
      if(score_table[o] > *max_score){
	*max_score = score_table[o];
	*max_row = row;
	*max_column = column;
      }
    }
  }
  

}
