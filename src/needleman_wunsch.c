#include <string.h>
#include <stdlib.h>
#include <string.h>
#include "needleman_wunsch.h"
#include <Rinternals.h>

#define INTRON_TYPE_N 3

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

const int left = 1;
const int up = 2;
const int diag = 3;

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
void init_ptr( int *ptr, int height, int width){  // const int left, const int up ){
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

  init_scores( score_table, height, width, gap_i, gap_e, tgaps_free );
  init_ptr( ptr_table, height, width); //, left, up );
  
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

      int left_score = score_table[o_left] + (((ptr_table[o_left] & ptr_mask) == left) ? gap_e : gap_i);
      int up_score = score_table[o_up] + (((ptr_table[o_up] & ptr_mask) == up) ? gap_e : gap_i);
      int diag_score = score_table[o_diag] + sub_table[ m_score_o ];

      // I used to make a vector of scores and use a function to calculate the max.
      // This is marginally faster.
      if(left_score > score_table[o]){
	score_table[o] = left_score;
	ptr_table[o] = left | (((( ptr_table[o_left] & left_mask ) >> left_shift) + 1) << left_shift );
	ptr_table[o] |=  (ptr_table[o_left] & up_mask);
      }
      if(up_score > score_table[o]){
	score_table[o] = up_score;
	ptr_table[o] = up | ( ((( ptr_table[o_up] & up_mask ) >> up_shift) + 1) << up_shift )  ;
	ptr_table[o] |=  ( ptr_table[o_up] & left_mask ); 
      }
      if(diag_score > score_table[o]){
	score_table[o] = diag_score;
	ptr_table[o] = diag | ( ((( ptr_table[o_diag] & left_mask ) >> left_shift) + 1) << left_shift );
	ptr_table[o] |= ( ((( ptr_table[o_diag] & up_mask ) >> up_shift) + 1) << up_shift ); 
      }
      if(score_table[o] > *max_score){
	*max_score = score_table[o];
	*max_row = row;
	*max_column = column;
      }
    }
  }
}


// This considers columns and rows different. It would be more appropriate to say:
// Extract non-overlapping alignments of sequence a (rows) to sequence b.
// Alignments can overlap in b (columns), but not in rows.
// In common with C++ iterators we take end to mean the position after the
// the last position; that is it is equivalent to the size of the table when
// starting counting at 0.
void extract_sw_alignments(const unsigned char *a, const unsigned char *b,
			   int *ptr_table, int *score_table, int height, int width,
			   int row_begin, int row_end, int col_begin, int col_end,
			   struct sw_alignment **sw_align, int min_width, int min_score)
{
  // make sure that the ranges are reasonble
  // hard coded for now; this ought to be changed.
  if( row_end - row_begin < min_width || col_end - col_begin < min_width)
    return;
  int max_score = 0;
  int max_row = 0;
  int max_column = 0;

  // default to there not being an alignment within the field:
  *sw_align = 0;
  if(row_end > height || col_end > width)
    return;
  
  for(int row=row_begin; row < row_end; ++row){
    for(int col=col_begin; col < col_end; ++col){
      int o = row + col * height;
      int l_moves = (ptr_table[o] & left_mask) >> left_shift;
      int u_moves = (ptr_table[o] & up_mask) >> up_shift;
      if( score_table[o] >= max_score && score_table[o] >= min_score &&
	  row - u_moves >= row_begin && col - l_moves >= col_begin ){
	max_score = score_table[o];
	max_row = row;
	max_column = col;
      }
    }
  }
  if(max_score == 0)
    return;
  // Then extract the number of positions in the residue. When doing this we can also do other things
  int al_length = 0;
  int row = max_row;
  int col = max_column;
  int o = row + col * height;
  int last_op = 0;
  int op_count = 0;
  while( (ptr_table[o] & ptr_mask) && (row > 0 || col > 0) ){
    if((ptr_table[o] & ptr_mask) != last_op){
      op_count++;
      last_op = (ptr_table[o] & ptr_mask);
    }
    col = (ptr_table[o] & left) ? col - 1 : col;
    row = (ptr_table[o] & up) ? row - 1 : row;
    o = row + col * height;
    al_length++;
  }
  // op_count lets us specify the cigar string; Note that the how we encode
  // the cigar string is a bit problematic. We could use a 14 bit represntation
  // giving us a maximum length of 16384, or a 30 bit representation, giving us
  // way more than we need. But this kind of string will not be very useful
  // in R; in any case it seems difficult to interrogate more than 31 bits
  // in R (which is probably related to the lack of unsigned integers).
  (*sw_align) = malloc( sizeof(struct sw_alignment) );
  struct sw_alignment *align = *sw_align;  // for ease of use
  memset( (void*)align, 0, sizeof(struct sw_alignment) );
  align->cigar_length = op_count;
  align->cigar_ops = malloc( sizeof(unsigned char) * op_count);
  align->cigar_n = malloc( sizeof(int) * op_count );
  memset( align->cigar_ops, 0, sizeof(unsigned char) * op_count );
  memset( align->cigar_n, 0, sizeof(int) * op_count );
  align->al_length = al_length;
  align->score = max_score;
  align->row_begin = row;
  align->row_end = max_row;
  align->col_begin = col;
  align->col_end = max_column;

  align->a_al = malloc( sizeof(char) * (al_length + 1) );
  align->b_al = malloc( sizeof(char) * (al_length + 1) );
  align->a_al[ al_length ] = 0;
  align->b_al[ al_length ] = 0;
    
  // and then we go through again and fill in the counts..
  int cigar_i = op_count; // decrement before use
  row = max_row;
  col = max_column;
  o = row + col * height;
  last_op = 0;
  int align_i = al_length-1;
  while( (ptr_table[o] & ptr_mask) && (row > 0 || col > 0) ){
    if((ptr_table[o] & ptr_mask) != last_op){
      cigar_i--;
      align->cigar_ops[cigar_i] = (ptr_table[o] & ptr_mask);
      last_op = (ptr_table[o] & ptr_mask);
    }
    align->a_al[ align_i ] = (ptr_table[o] & up) ? a[ row-1 ] : '-';
    align->b_al[ align_i ] = (ptr_table[o] & left) ? b[col-1] : '-';
    align_i--;
    row = (ptr_table[o] & up) ? row - 1 : row;
    col = (ptr_table[o] & left) ? col - 1 : col;
    o = row + col * height;
    align->cigar_n[cigar_i]++;
  }
  // and then we recurse to the four remaining quadrants.. 
  extract_sw_alignments(a, b, ptr_table, score_table, height, width,
			row_begin, align->row_begin, col_begin, col_end,
			&(align->top), min_width, min_score);
  extract_sw_alignments(a, b, ptr_table, score_table, height, width,
			align->row_end+1, row_end, col_begin, col_end,
			&(align->bottom), min_width, min_score);
}

void free_sw_alignments( struct sw_alignment *align )
{
  if(!align)
    return;
  free_sw_alignments( align->top );
  free_sw_alignments( align->bottom );
  free( align->cigar_ops );
  free( align->cigar_n );
  free( align->a_al );
  free( align->b_al );
  free( align );
}

void count_sw_alignments( struct sw_alignment *align, int *n)
{
  if(!align)
    return;
  (*n)++;
  count_sw_alignments( align->top, n );
  count_sw_alignments( align->bottom, n );
}

void harvest_sw_aligns( struct sw_alignment *align, int *align_table,
			int **cigar_ops, int **cigar_n, int *cigar_lengths,
			int *i, int aligns_n, char **a_al, char **b_al ){
  if(!align || *i >= aligns_n)
    return;

  // copy the data to the i th position;
  align_table[ (*i) + aligns_n * 0] = align->row_begin;
  align_table[ (*i) + aligns_n * 1] = align->row_end;
  align_table[ (*i) + aligns_n * 2] = align->col_begin;
  align_table[ (*i) + aligns_n * 3] = align->col_end;
  align_table[ (*i) + aligns_n * 4] = align->score;
  align_table[ (*i) + aligns_n * 5] = align->al_length;
  cigar_lengths[*i] = align->cigar_length;
  cigar_n[*i] = malloc( sizeof(int) * align->cigar_length );
  cigar_ops[*i] = malloc( sizeof(int) * align->cigar_length );
  for(int j=0; j < align->cigar_length; ++j){
    cigar_n[*i][j] = align->cigar_n[j];
    cigar_ops[*i][j] = align->cigar_ops[j];
  }
  a_al[(*i)] = align->a_al;
  b_al[(*i)] = align->b_al;
  // increment the operator..
  (*i)++;
  harvest_sw_aligns( align->top, align_table, cigar_ops, cigar_n, cigar_lengths, i, aligns_n, a_al, b_al );
  harvest_sw_aligns( align->bottom, align_table, cigar_ops, cigar_n, cigar_lengths, i, aligns_n, a_al, b_al );
}

// align_transcript
// Takes an intron-marked transcript sequences and aligns to a genomic loci.
// Introns are marked by a special character (I), which has a high mismatch penalty
// but no gap insertion (in the genomic sequence).
// Since genomic loci can be long we must implement a memory efficient form of the
// the alignment which only returns the aligned regions; i.e. the presumptive exons.

struct transcript_alignment {
  size_t intron_n;
  int *intron_beg, *intron_end;
  size_t exon_n; // should be intron_n + 1
  int *exon_beg, *exon_end;
};


// This function is very complex and should almost certainly be rewritten.. 
// Using more sane data structures to reduce the very large number of malloc() calls.
// 
// THIS function needs to be broken down into more reasonable structs and functions
// at the moment it is far too complicated and has far too small a chance of actually
// working. But I need to do this at a time when I can give more time to the problem.
//
// intron_char is the character used to represent introns in tr_seq
// sub_table is the substitution table, al_offset is the alphabet offset
// and al_size is the alphabet size
// gap_i and gap_e are penalties for inserting and extending a gap
// intron_loss is the penalty for a loss of an intron in the genome sequence.
/* #if FALSE */
/* void align_transcript( struct transcript_alignment *tr_alignment, */
/* 		       const char *tr_seq, const char *g_seq, */
/* 		       int gap_i, int gap_e, int intron_loss, char intron_char, */
/* 		       int *sub_table, int al_offset, int al_size ) */
/* { */
/*   // We consider R-based matrix coordinates (column major) and that the tr_seq is vertical */
/*   // with the g_seq along the top. */
/*   // To find the coordinates of the highest scoring alignment we need only to keep track of the previous */
/*   // set of scores and alignment types; We will remember only the highest scoring cell for any given column */
/*   // along with the sum of the operations required to get to that cell. */
/*   // We also need to remember the starts and lenths of introns; For each of these maxima. But we only need */
/*   // to remember that for the highest scoring position. */
  
/*   // First determine the number of introns; */
/*   int intron_n = 0; */
/*   size_t tr_seq_l = 0; */
/*   for( char nt=tr_seq; *nt > 0; nt++ ){ */
/*     ++tr_seq_l; */
/*     if( *nt == intron_char ) */
/*       ++intron_n; */
/*   } */
/*   int exon_n = intron_n + 1; */
/*   // these are the scores and the states associated with the previous column */
/*   // Since at the end we need to remember the beginnings and lengths of all */
/*   // potential introns we make a table for all of them */
/*   int *p_scores = malloc( sizeof(int) * (1 + tr_seq_l) ); */
/*   char *p_ops = malloc( sizeof(char) * (1 + tr_seq_l) ); // see below for encoding */
/*   // how many insertions and deletions (referring to the transcript) */
/*   // for a given score; */
/*   // We need one per intron..  */
/*   // and we probably need to keep track of the number of matches */
/*   // as well..  */
/*   int *insert_n = malloc( sizeof(int) * tr_seq_l * intron_n); */
/*   int *delete_n = malloc( sizeof(int) * tr_seq_l * intron_n); */
/*   int *match_n =  malloc( sizeof(int) * tr_seq_l * intron_n); */

/*   // information about intron positions; also refer to the previous column */
/*   // of scores */
/*   int *intron_pos = malloc( sizeof(int) * tr_seq_l * intron_n ); */
/*   int *intron_length = malloc( sizeof(int) * tr_seq_l * intron_n ); */
/*   int *intron_i = malloc( sizeof(int) * tr_seq_l ); */

/*   // these should all be set to 0 to start with. */
/*   memset( (void*)p_scores, 0, sizeof(int) * (1 + tr_seq_l) ); */
/*   memset( (void*)p_operations, 0, sizeof(char) * (1 + t_seq_l) ); */
/*   memset( (void*)insert_n, 0, sizeof(int) * tr_seq_l * intron_n ); */
/*   memset( (void*)delete_n, 0, sizeof(int) * tr_seq_l * intron_n ); */
/*   memset( (void*)match_n, 0, sizeof(int) * tr_seq_l * intron_n ); */
/*   memset( (void*)intron_pos, 0, sizeof(int) * tr_seq_l * intron_n ); */
/*   memset( (void*)intron_length, 0, sizeof(int) * tr_seq_l * intron_n ); */
/*   memset( (void*)intron_i, 0, sizeof(int) * tr_seq_l ); */

/*   // The maximum score along the genomic locus */
/*   int max_score; */
/*   int max_i;     // (the index or position of the score) */
/*   int *max_insert_n = malloc(sizeof(int) * intron_n); */
/*   int *max_delete_n = malloc(sizeof(int) * intron_n); */
/*   int *max_intron_pos = malloc(sizeof(int) * intron_n); */
/*   int *max_intron_length = malloc(sizeof(int) * intron_n); */
/*   int max_intron_i; */
  
/*   size_t g_seq_l = strlen( g_seq ); */
/*   for(int i=0; i < g_seq_l; ++i){ */
/*     int pp_score = 0; */
/*     int pp_op = 0; */
/*     for(int j=1; j <= tr_seq_l; ++j){ */
/*       // Determine the score for the virtual cell in row j and column (i+1) */
/*       // based on the scores in insert_n; */
/*       int cell_max = 0; */
/*       int cell_op = 0; */
/*       // operations are:  */
/*       // insert (move down with gap) */
/*       // intron_insertion (move down, tr seqeunce is I), heavy penalty */
/*       // delete (move right with gap) */
/*       // intron (move right, tr sequence is I), no penalty */
/*       int ins_penalty = (tr_seq[j-1] == intron_char) ? intron_loss : */
/* 	pp_ops == 1 ? gap_e : gap_i; */
/*       int del_penalty = (tr_seq[j-1] == intron_char) ? 0 : */
/* 	p_ops[j] == 2 ? gap_e : gap_i; */
/*       int match_penalty = sub_table[ (g_seq[i]-al_offset) + (tr_seq[j] - al_offset) * al_size ]; */
/*       if(p_scores[j-1] + ins_penalty > cell_max ){ */
/* 	cell_max = pp_score + ins_penalty; */
/* 	cell_op = 1; */
/*       } */
/*       if(p_scores[j] + del_penalty > cell_max ){ */
/* 	cell_max = p_scores[j] + del_penalty; */
/* 	cell_op = 2; */
/*       } */
/*       if(p_scores[j-1] + match_penalty > cell_max ){ */
/* 	cell_max = p_scores[j-1] + del_penalty; */
/* 	cell_op = 4; */
/*       } */
/*       if( (tr_seq[j-1] == intron_char) ) */
/* 	cell_op |= 8; */
/*       p_scores[j] = cell_max; */
/*       p_ops[j] = cell_op; */
/*       pp_score = cell_max; */
/*       pp_op = cell_op; */
/*       // then update tha maximum and and all the counts; */
/*       if(cell_score == 0){  // no alignment cannot be maximal..  */
/* 	intron_i[j] = 0; */
/* 	// the following is very wasteful and should not be necessary most of the */
/* 	// time */
/*       } */

/*     } */
/*   } */
  

/* } */
/* #endif */

// The following code replaces the transcript_align function above
// through restructing the data into logical units..

struct cig_op *init_cig_op(unsigned char op, struct cig_op *parent){
  struct cig_op *node = malloc( sizeof(struct cig_op) );
  node->parent = parent;
  node->data = op; // ops are 1: left, 2 up, 3 diagonal
  if(parent && op)
    parent->data |= (1 << (1 + op));
  return(node);
}

unsigned char c_op(struct cig_op *co){
  return(co->data & 3);
}

void clear_ops(struct cig_op *node, unsigned char op){
  if(!node)
    return;
  unsigned char mask = (1 << (1+op));
  if(op)
    node->data &= (~mask);
  // if any of bits 3-5 are set.. 
  if(node->data & 0x1C) // some other cell pointing at this one
    return;
  if(node->parent)
    clear_ops(node->parent, node->data & 3);
  free(node);
}


struct cig_data *init_cig_data(unsigned int l){
  struct cig_data *cd = malloc(sizeof(struct cig_data) * l);
  memset( cd, 0, sizeof(struct cig_data) * l);
  return(cd);
}

// clear the values in cd;
// if no children, then also clear the cig_op linked list
// if cd->max_score is better than best->max_score 
// AND (op OR no children)
//     copy cd -> best
//     clear cd to 0
//     do not clear the cig_ops
//
void clear_cig_data(struct cig_data *cd, struct cig_data *best){
  if(cd->ops && (cd->ops->data >> 2 == 0) && cd->max_score <= best->max_score)
    clear_ops(cd->ops, 0);
  if(cd->max_score > best->max_score && (cd->ops->data >> 2 == 0)){
    clear_ops(best->ops, 0);
    memcpy(best, cd, sizeof(struct cig_data));
    // and DO NOT clear cd->ops!
  }
  memset(cd, 0, sizeof(struct cig_data)); // this could be changed to a single memset op at the end of the inner loop
}

int extract_swi_alignment( struct cig_data *cd, const unsigned char *tr, const unsigned char *gen,
					   struct sw_alignment *alignment){
  // first find the top scoring position and then count from that:
  // 1. How many distinct cigar operations are encoded
  // 2. How many distinct moves from the top sc
  struct cig_op *max_op = 0;
  struct cig_op *op = cd->ops;
  int pos_t = cd->pos_t;
  int pos_g = cd->pos_g;
  while(op && (pos_t > cd->max_t || pos_g > cd->max_g)){
    // decrement pos_t & pos_g appropriately here
    pos_t -= (op->data & 2) ? 1 : 0;
    pos_g -= (op->data & 1) ? 1 : 0;
    op = op->parent;
  }
  // if not op at this point then something has gone decidedly wrong
  if(!op){
    Rprintf("Major error, extract swi_alignment obtained NULL op at alignment max\n");
    return(1); 
  }
  if(pos_t != cd->max_t || pos_g != cd->max_g){
    Rprintf("pos_t or pos_g has bad value\n");
    return(2);
  }
  int op_n = 0;  // all operations
  int op_dn = 0; // distinct operations
  unsigned char last_op = 0;
  max_op = op;
  while(op){
    op_n++;
    op_dn += (last_op != (op->data & 3)) ? 1 : 0;
    last_op = op->data & 3;
    op = op->parent;
  }
  // then we can set up the structures that we need to hold the alignment;
  alignment->score = cd->max_score;
  // row and column are one plus the position; 
  alignment->row_begin = cd->beg_t + 1;
  alignment->row_end = cd->max_t + 1;
  alignment->col_begin = cd->beg_g + 1;
  alignment->col_end = cd->max_g + 1;
  alignment->al_length = op_n;
  alignment->a_al = malloc( op_n + 1 );
  alignment->b_al = malloc( op_n + 1 );
  alignment->a_al[op_n] = 0;
  alignment->b_al[op_n] = 0;
  alignment->cigar_length = op_dn;
  alignment->cigar_n = malloc(op_dn * sizeof(int));
  memset( alignment->cigar_n, 0, sizeof(int) * op_dn);
  alignment->cigar_ops = malloc(op_dn);
  alignment->top = 0;
  alignment->bottom = 0;
  // then simply go through and set values appropriately
  op = max_op;
  last_op = 0;
  while(op){
    op_n--;
    // the following is not the most efficient way as
    // we end up with more conditionals than necessary
    // but it seems easier to understand
    alignment->a_al[op_n] = (op->data & 2) ? tr[pos_t] : '-';
    alignment->b_al[op_n] = (op->data & 1) ? gen[pos_g] : '-';
    pos_t -= (op->data & 2) ? 1 : 0;
    pos_g -= (op->data & 1) ? 1 : 0;
    if(last_op != (op->data & 3)){
      op_dn--;
      alignment->cigar_ops[op_dn] = op->data & 3;
    }
    alignment->cigar_n[op_dn]++;
    last_op = (op->data & 3);
    op = op->parent;
  }
  return(0);
}

// WARNING
// This function may not work, or if it does, it is very slow
// The function attemptes to keep cigar representations of alignments
// as linked lists (to avoid keeping a table of the pointers for the full
// alignment matrix). Not surprisingly, this turns out to be really computationally
// expensive, and not really worth it.
// It might be possible to make this work, through using preallocation and a specialised
// memory mapper, but that is probably not worth it.
// Instead, use the function further below: align_transcript_bitp()
// Which simply uses two bits per cell and keeps the whole table. Not very elegant, but
// works much faster.
struct sw_alignment *align_transcript_to_genome( const unsigned char *tr, const unsigned char *gen,
						int tr_l, int gen_l, int gap_i, int gap_e, int intron_p,
						int *sub_table, int al_offset, int al_size, char intron_char){
  int al_buffer_size = 32;
  int left_max_gap = -1 * (1 + abs(intron_p - gap_i) / gap_e);
  // use to allow long gaps upwards as well... necessary if an alignment has sequence that doesn't fit..
  // left gap insertion penalties are 0 after left_max_gap operations.
  // to introns. 
  // We need to have a vector of cig_data the same size as the transcript. We will
  // trust that we have the correct size of these.
  // this represents the current alignments, so we will call it aligns
  struct cig_data *left = init_cig_data( (unsigned int)tr_l + 1 );
  struct cig_data *right = init_cig_data( (unsigned int)tr_l + 1 );
  struct cig_data best_alignment;
  memset(&best_alignment, 0, sizeof(struct cig_data));
  // everything starts of as zero
  for(int g_i=0; g_i < gen_l; ++g_i){
    for(int tr_i=0; tr_i < tr_l; ++tr_i){
      int row = tr_i + 1;
      // The offset that gives us the match / mismatch penalty
      int m_score_o = (tr[tr_i] - al_offset) + (gen[g_i] - al_offset) * al_size;
      unsigned int left_op = left[row].ops ? left[row].ops->data & 3 : 0;
      unsigned int up_op = right[row-1].ops ? right[row-1].ops->data & 3 : 0;
      // left scores never decrease below (score - intron_p)
      // the lengt of a left gap affects the diagonal score rejoining an alignment
      // if we have more than left.n consecutive gap positions then we do not penalise further extension.
      int left_score = left[row].score + ((left_op & 3) == 1 ? gap_e : (left[row].left_n > left_max_gap ? 0 : gap_i));
      int up_score = right[row-1].score + ((up_op & 3) == 2 ? gap_e : (right[row-1].up_n > left_max_gap ? 0 : gap_i));
      // we might want to consider modifying the diagonal score if the previous character of the
      // transcript was an intron character; possibly to penalise introns that are too short! We could also
      // set a penalty here based on the intron end characters (eg. GT-AC, AT-AC, and so on to avoid
      // incorrectly spliced ones... 
      int diag_score = left[row-1].score + sub_table[ m_score_o ];
      // we might want to make left_score and diag_score have no penaties at all if we are on an intron row
      // although, arguably, maybe the up_score should be the intron_penalty as not having an intron here
      // is the same as having lost one..
      //      Rprintf("%d  %d  %d  %d\n", best_alignment.score, left_score, up_score, diag_score);
      if(tr[tr_i] == intron_char){
	up_score = right[row-1].score + intron_p;  // new intron
	left_score = left[row].score;  // no penalty at all.. this is the novelty in this alignment.
	diag_score = left[row-1].score; // no penalty at all. We have to join in some way.. 
      }
      // We could also consider to modify the diagonal score so that rejoining from an intron will take
      // an indel based penalty if the gap is too small for the intron.
      // The default value is 0, with no cigar...
      int scores[4] = {0, left_score, up_score, diag_score};
      struct cig_data *parents[4] = {0, left+row, right+row-1, left+row-1};
      int max_i = 0;
      struct cig_op *parent = 0;
      for(int i=1; i < 4; ++i)
	max_i = (scores[i] > scores[max_i]) ? i : max_i;
      // if max_i is 0, then we don't actually need to do very much as it has no information
      //      Rprintf("tr_i: %d tr_l: %d row: %d  gen: %d\n ", tr_i, tr_l, row, g_i);
      if(max_i == 0){
	memset( right + row, 0, sizeof(struct cig_data) );
      }else{
	// most of the fields should default to the parent one
	memcpy(right + row, parents[max_i], sizeof(struct cig_data));
	right[row].ops = init_cig_op(max_i, parents[max_i]->ops);
	right[row].score = scores[max_i];
	// if the first one; that is the parent is empty
	if(parents[max_i]->ops == 0){
	  right[row].beg_t = tr_i;
	  right[row].beg_g = g_i;
	}
	if(right[row].score > right[row].max_score){
	  right[row].max_score = right[row].score;
	  right[row].max_t = tr_i;
	  right[row].max_g = g_i;
	}
	if(max_i == 1)
	  right[row].left_n = parents[max_i]->left_n + 1;
	else
	  right[row].left_n = 0;
	if(max_i == 2)
	  right[row].up_n = parents[max_i]->up_n + 1;
	else
	  right[row].up_n = 0;
      }
      // we always set the positions of the cig_data
      right[row].pos_t = tr_i;
      right[row].pos_g = g_i;
      // Clear the data from above left if no children
      clear_cig_data(&left[row-1], &best_alignment);
    }
    //    Rprintf("Got here, clearing tr_l\n");
    clear_cig_data(&left[tr_l], &best_alignment);
    // swap left and right
    struct cig_data *tmp = left;
    left = right;
    right = tmp;
  }
  // and go through the last column, clearing every element. Note at this point the
  // last column is left as they are swapped at the end. 
  for(int row=1; row <= tr_l; ++row)
    clear_cig_data( &left[row], &best_alignment );
  // Free up left and right:
  free(left);
  free(right);
  // Now the best alignment should be present in the best_alignment and we can extract that
  struct sw_alignment *alignment = malloc(sizeof(struct sw_alignment));
  memset(alignment, 0, sizeof(struct sw_alignment));
  int error = extract_swi_alignment(&best_alignment, tr, gen, alignment);
  clear_ops( best_alignment.ops, 0 );
  return(alignment);
}

struct bit_ptr init_bit_ptr(unsigned int nrow, unsigned int ncol){
  struct bit_ptr ptr;
  ptr.nrow=nrow;
  ptr.ncol=ncol;
  unsigned int ptr_n = nrow * ncol;
  unsigned int n_bytes = (ptr_n % 4) ? 1 + ptr_n / 4 : ptr_n / 4;
  ptr.data = malloc( ptr_n );
  memset(ptr.data, 0, ptr_n );
  return(ptr);
}

void free_bit_ptr(struct bit_ptr *ptr){
  free( ptr->data );
}

// seems like this would be a good time to define a macro or inline function
unsigned char get_ptr(struct bit_ptr *ptr, unsigned int row, unsigned int col){
  size_t o = row + col * ptr->nrow;
  return( 3 & ptr->data[ o / 4 ] >> 2 * (o % 4) );
}

void set_ptr(struct bit_ptr *ptr, unsigned int row, unsigned int col, unsigned char val){
  // the simple expression is:
  size_t o = row + col * ptr->nrow;
  ptr->data[ o / 4 ] |= (val << 2 * (o % 4) );
  // but this only works if the original ptr data is 0, otherwise it's not a setting operation.
}

int extract_bitp_alignment( const unsigned char *tr, const unsigned char *gen, 
			    struct bit_ptr *ptr, int max_row, int max_col, int score,
			    struct sw_alignment *al){
  // First count the number of operations needed to store the data;
  int op_n = 0; // total number of operations
  int dop_n = 0; // number of distinct operations
  unsigned char op = 0;
  unsigned char last_op = 0;
  int row = max_row;
  int col = max_col;
  while( (op = get_ptr(ptr, row, col)) ){
    op_n++;
    if(op != last_op)
      dop_n++;
    if(op & 1)
      col--;
    if(op & 2)
      row--;
    last_op = op;
  }
  al->a_al = malloc(op_n + 1);
  al->b_al = malloc(op_n + 2);
  al->a_al[op_n] = 0;
  al->b_al[op_n] = 0;
  al->score = score;
  al->row_begin = row+1;
  al->row_end = max_row;
  al->col_begin = col+1;
  al->col_end = max_col;
  al->al_length = op_n;
  al->cigar_length = dop_n;
  al->cigar_ops = malloc(dop_n);
  al->cigar_n = malloc(sizeof(int) * dop_n);
  memset( al->cigar_n, 0, sizeof(int) * dop_n);
  al->top=0;
  al->bottom=0;
  // and then set the cigar and things appropriatly
  row = max_row;
  col = max_col;
  last_op = 0;
  op = 0;
  // I should possibly be checking for a corrupted pointer table here
  // but... 
  while( (op = get_ptr(ptr, row, col)) && op_n > 0 ){
    op_n--;
    al->a_al[op_n] = (op & 2) ? tr[row-1] : '-';
    al->b_al[op_n] = (op & 1) ? gen[col-1] : '-';
    if(op != last_op){
      dop_n--;
      al->cigar_ops[dop_n] = op;
    }
    al->cigar_n[dop_n]++;
    if(op & 1)
      col--;
    if(op & 2)
      row--;
    last_op = op;
  }
  return(0);
}

struct sw_alignment *align_transcript_bitp(const unsigned char *tr, const unsigned char *gen,
					   int tr_l, int gen_l, int gap_i, int gap_e, int intron_p,
					   int *sub_table, int al_offset, int al_size, char intron_char)
{
  int intron_gap_max = -1 * (1 + abs(intron_p - gap_i) / gap_e);
  int *left_scores = malloc(sizeof(int) * (tr_l + 1));
  int *right_scores = malloc(sizeof(int) * (tr_l + 1));
  memset(left_scores, 0, sizeof(int) * (tr_l + 1));
  memset(right_scores, 0, sizeof(int) * (tr_l + 1));
  // I may at some point make the left_scores and right_scores into data structures that also
  // contain the starting point of the alignment and the counts of current left or up ptrs
  // to modify the scoring;
  int *left_count = malloc(sizeof(int) * (tr_l + 1));
  int *up_count = malloc(sizeof(int) * (tr_l + 1));
  memset(left_count, 0, sizeof(int) * (tr_l + 1));
  memset(up_count, 0, sizeof(int) * (tr_l + 1));
  // and the pointers.. 
  struct bit_ptr ptr = init_bit_ptr( tr_l + 1, gen_l + 1);
  // and the maximum position;
  int max_score = 0;
  int max_col = 0;
  int max_row = 0;
  // This should be enough...
  for(int g_i=0; g_i < gen_l; ++g_i){
    int col = g_i + 1;
    for(int tr_i=0; tr_i < tr_l; ++tr_i){
      int row = tr_i + 1;
      // then set right_scores[i]
      // and update left_count and up_count to indicate the previous count..
      int m_score_o = (tr[tr_i] - al_offset) + (gen[g_i] - al_offset) * al_size;
      // I might not need these pointers since I have the counts, let me think..
      
      int left_score = left_scores[row] + ((left_count[row] == 0) ? gap_i : (left_count[row] < intron_gap_max ? gap_e : 0));
      int up_score = right_scores[row-1] + ((up_count[row-1] == 0) ? gap_i : (up_count[row-1] < intron_gap_max ? gap_e : 0));
      int diag_score = left_scores[row-1] + sub_table[ m_score_o ];

      // intron character only allowed in tr
      if(tr[tr_i] == intron_char){
	up_score = right_scores[row-1] + intron_p;  // new intron
	left_score = left_scores[row];  // no penalty at all.. this is the novelty in this alignment.
	diag_score = left_scores[row-1]; // no penalty at all. We have to join in some way.. 
      }
      int scores[4] = {0, left_score, up_score, diag_score};
      int max_i = 0;
      for(int i=1; i < 4; ++i)
	max_i = (scores[i] > scores[max_i]) ? i : max_i;
      // and then simply:
      left_count[row] = (max_i == 1) ? left_count[row] + 1 : 0;
      // up_count[row] + 1 --> up_count[row-1] + 1
      // 19-06-23 lmj. Looks like previous version was a bug.
      up_count[row] = (max_i == 2) ? up_count[row-1] + 1 : 0;
      right_scores[row] = scores[max_i];
      set_ptr( &ptr, row, col, (unsigned char)max_i );
      // update the maximum if appropriate;
      if(right_scores[row] > max_score){
	max_score = right_scores[row];
	max_col = col;
	max_row = row;
      }
    }
    int *tmp = left_scores;
    left_scores = right_scores;
    right_scores = tmp;
  }
  struct sw_alignment *alignment = malloc(sizeof(struct sw_alignment));
  memset(alignment, 0, sizeof(struct sw_alignment));
  int error = extract_bitp_alignment( tr, gen, &ptr, max_row, max_col, max_score, alignment);
  // clear all the allocated memory:
  free(left_scores);
  free(right_scores);
  free(left_count);
  free(up_count);
  free_bit_ptr(&ptr);
  return(alignment);
}


/* // returns the nth class */
/* int check_intron_entry(const char *seq, const char **allowed, unsigned char n){ */
/*   for(unsigned char i=0; i < n; ++i){ */
/*     if(strncmp(seq, allowed[i], 2) == 0) */
/*       return(i+1); */
/*   } */
/*   return(0); */
/* } */

/* // returns the nth class if intron merge is OK */
/* // This may not be necessary if we remember which intron class we have entered. */
/* // But for that we will need another set of variables. Unsigned char ?, with one for each */
/* // entry...  */
/* int check_intron_close(const char *seq1, const char *seq2, const char **allowed_start, const char **allowed_end, unsigned char n){ */
/*   for(unsigned char i=0; i < n; ++i){ */
/*     if(strncmp(seq1, allowed_start[i], 2) == 0 && strncmp(seq2, allowed_end[i]) == 0) */
/*       return(i+1); */
/*   } */
/*   return(0); */
/* } */

int check_intron(const char *seq1, const char *seq2, const char **allowed_start, const char **allowed_end, unsigned char n){
  for(unsigned char i=0; i < n; ++i){
    if(strncmp(seq1, allowed_start[i], 2) == 0 && strncmp(seq2, allowed_end[i], 2) == 0)
      return(i+1);
  }
  return(0);
}

// Contains two columns of scores, plus two columns of gap counts
// for the current (right) and previous (left) columns of the positions
// The vertical gap counts (sequence in transcript not in 
struct al_info {
  // pointers to the actual sequences
  // and their lenghts.
  const unsigned char *tr;
  const unsigned char *gen;
  unsigned int tr_l, gen_l;
  // the character associated with an intron
  char intron_char;
  // The two column of scores; these get swapped at the
  // end of an assignment
  float *left_scores;
  float *right_scores;
  // For each of these vectors we also keep track of the number of the
  // size of any gaps to the corresponding position.
  // tr_gap corresponds to left arrows and are presumptive introns.
  unsigned int *left_tr_gap, *right_tr_gap;
  // gen_gap corresponds to up arrows and represent sequence in tr that
  // is not aligned to the genome (insertions in the transcript).
  // These are not expected but may be present if part of an exon is lost
  // in the species of the genome sequence.
  unsigned int *left_gen_gap, *right_gen_gap;
  // gaps in genome sequence (up arrows) can cross intron positions in the transcript
  // If they do we need to count them in order to apply the appropriate penalties. 
  unsigned int *left_int_x, *right_int_x;
};

struct al_info init_al_info(const unsigned char *tr, const unsigned char *gen,
			    unsigned int tr_l, unsigned int gen_l, char intron_char){
  struct al_info al;
  al.tr = tr;
  al.gen = gen;
  al.tr_l = tr_l;
  al.gen_l = gen_l;
  al.intron_char = intron_char;
  // calloc clears the allocated memory setting it to 0.
  // Note that this could be done with a single allocation
  // of a block of memory from which we could take a subset.
  al.left_scores = calloc( tr_l + 1, sizeof(float) );
  al.right_scores = calloc( tr_l + 1, sizeof(float) );
  al.left_tr_gap =  calloc( tr_l + 1, sizeof(unsigned int) );
  al.right_tr_gap =  calloc( tr_l + 1, sizeof(unsigned int) );
  al.left_gen_gap =  calloc( tr_l + 1, sizeof(unsigned int) );
  al.right_gen_gap =  calloc( tr_l + 1, sizeof(unsigned int) );
  // and the intron crosses for gen_gaps
  // Note that these really _dont_ need to be 32 bit ints
  al.left_int_x = calloc( tr_l + 1, sizeof(unsigned int) );
  al.right_int_x = calloc( tr_l + 1, sizeof(unsigned int) );
  return(al);
}


// It strikes me that we do not actually need to swap; 
void swap_left_right( struct al_info *al ){
  float *tmp = al->left_scores;
  al->left_scores = al->right_scores;
  al->right_scores = tmp;
  unsigned int *utmp;
  // tr gaps
  utmp = al->left_tr_gap;
  al->left_tr_gap = al->right_tr_gap;
  al->right_tr_gap = utmp;
  // gen gaps
  utmp = al->left_gen_gap;
  al->left_gen_gap = al->right_gen_gap;
  al->right_gen_gap = utmp;
  // int_x counts
  utmp = al->left_int_x;
  al->left_int_x = al->right_int_x;
  al->right_int_x = utmp;
}

void free_al_info(struct al_info *al){
  free( al->left_scores );
  free( al->right_scores );
  free( al->left_tr_gap );
  free( al->right_tr_gap );
  free( al->left_gen_gap );
  free( al->right_gen_gap );
  free( al->left_int_x );
  free( al->right_int_x );
}

// A struct to hold the set of penalties that we wish to use:
// defined in header file.
/* struct al_penalty { */
/*   int *sub_table; */
/*   int al_offset; */
/*   int al_size; */
/*   char intron_char; */
/*   float gap_i; */
/*   // gap_e_mod modifies the extension: the bigger the cheaper the extension as it divides */
/*   float gap_e_mod;  */
/*   float fshift[3]; */
/*   float intron_new; */
/*   float intron_bad; */
/*   int min_intron_size; */
/*   // intron starts and ends: */
/*   int intron_type_n; */
/*   // precalculated cumulative gap values for 0 -> min_intron_size; */
/*   // (size of array = min_intron_size + 1); */
/*   float *gap; */
/*   const char **intron_start; */
/*   const char **intron_end; */
/* }; */

// this is just a holder; it should not allocate or free anything
struct al_penalty init_al_penalty(int *sub_table, int al_offset, int al_size, char intron_char,
				  float gap_i, float gap_e_mod, float fshift, float intron_new, float intron_bad, int min_intron_size,
				  int intron_type_n, const char **intron_start, const char **intron_end
				  ){
  struct al_penalty pen;
  pen.sub_table = sub_table;
  pen.al_offset = al_offset;
  pen.al_size = al_size;
  pen.intron_char = intron_char;
  pen.gap_i = gap_i;
  pen.gap_e_mod = gap_e_mod;
  pen.fshift[0] = -fshift;
  pen.fshift[1] = fshift;
  pen.fshift[2] = 0;
  pen.intron_new = intron_new;
  pen.intron_bad = intron_bad;
  pen.min_intron_size = min_intron_size;
  pen.intron_type_n = intron_type_n;
  pen.intron_start = intron_start;
  pen.intron_end = intron_end;
  pen.gap = 0;
  if(min_intron_size > 0){
    pen.gap = malloc( sizeof(float) * (1 + min_intron_size));
    pen.gap[0] = 0;
    for(unsigned int i=1; i <= (unsigned int)min_intron_size; ++i)
      pen.gap[i] = pen.gap[i-1] + gap_penalty(i, (float)min_intron_size, gap_i, gap_e_mod, pen.fshift);
  }
  return(pen);
}

void free_al_penalty(struct al_penalty *pen){
  free(pen->gap);
}
 
// Individual functions to calculate the scores: left, up, diagonal
// and gap-closing for horizontal gaps (potential introns)

// gap_l can be either al->right_gen_gap or al->left_gen_gap depending on which cell is tested
// row is the row from which the closing is done. Hence we check row-1, rather than row.
// That is probably a bad choice as it conflicts with vertical_gap_close_penalty where the function
// can't assume row-1 and the caller has to be explicit.
static inline float horizontal_gap_close_penalty(struct al_info *al, struct al_penalty *pen, int row, unsigned int *gap_lengths, int g_i){
  if(row < 2 || gap_lengths[row-1] < 6) 
    return(0);
  // we will make use of a bitwise flag with 3 bits defining the following:
  // bit 1 : gap_l >= pen->min_intron_size
  // bit 2 : preceding row was an intron
  // bit 3 : intron_class > 0 (i.e. looks like a good intron)
  unsigned int flag = 0;
  unsigned int gap_l = gap_lengths[row-1];
  flag |= (al->tr[row-2] == pen->intron_char);
  flag |= ((gap_l >= pen->min_intron_size) << 1);
  flag |= (check_intron(al->gen + g_i - gap_l, al->gen + g_i - 1, pen->intron_start, pen->intron_end, pen->intron_type_n) << 2);
  // The flag represents the following table (Most significant bit first)
  // intron_class   intron_size    intron_row     penalty             reason
  //            0             0             0     0                  counts as normal gap: already calculatd
  //            0             0             1     gap                current 0 penalty, but should count as normal gap
  //            0             1             0     0                  counts as normal gap of intron_min_size; already calculated
  //            0             1             1     intron_bad / gap   current cost 0, but doesn't look like intron
  //            1             0             0     0                  doesn't look like intron, gap cost already included
  //            1             0             1     intron_bad / gap   or some kind of function depending on gap_l
  //            1             1             0     intron_new         looks like new intron
  //            1             1             1     0                  looks like good intron
  // Note that in some of these case we may need to return a negative number. 110 (intron class, size) should
  // have a total penalty of intron_new, but will already have a gap penalty that may be larger than the
  // than the intron_new penalty. And 001, or 011, should have a value for the gap length; that suggests
  // that we should pre-calculate a gaps vector up to intron_min_size that we can use here.
  float penalty = 0;
  switch(flag){
  case 0:  // 0 0 0
    penalty = 0;
    break;
  case 1:  // 0 0 1
    // intron row, but too short, 
    penalty = 0;
    //    penalty = pen->gap[ gap_l ];
    break;
  case 2:  // 0 1 0 intron_size only
    penalty = 0;
    break;
  case 3:  // 0 1 1 intron row and intron_size, but bad beg / end
    // penalty = 0;
    penalty = pen->intron_bad;
    //    penalty = pen->gap[ pen->min_intron_size ];
    break;
  case 4:  // 1 0 0 intron class only
    penalty = 0;
    break;
  case 5:  // 1 0 1
    // intron class, intron row, but too short? 
    // penalty = 0;
    penalty = pen->gap[ gap_l ];
    break;
  case 6:  // 1 1 0
    // intron class, intron size, but not intron row
    // provide a positive penalty (the others are 0 or negative)
    penalty = 0;
    //    penalty = pen->intron_new - pen->gap[ pen->min_intron_size ];
    break;
  case 7: // 1 1 1
    // Everything looks like an intron !
    penalty = 0;
    break;
  default:
    penalty = 0; // should not be possible!
  }
  return(penalty);
}

// gap_l is the gap at this cell; not the previous one that is recorded!
// fshift should be positive and should reduce the penalty as it should not accumulate
// over the two first 
// complains if I make this inline.. 
static inline float gap_penalty(unsigned int gap_l, float min_is, float gap_i, float mod, float fshift[3]){
  return( (gap_l == 1 ? gap_i : 0) + logf( (float)gap_l / min_is )/mod + fshift[gap_l % 3] );
}

float score_from_left(struct al_info *al, struct al_penalty *pen, int row){
  // If we are on an intron row, then there is no penalty:
  if( al->tr[row-1] == pen->intron_char || al->left_tr_gap[row] >= pen->min_intron_size)
    return(al->left_scores[row]);
  // hardcode these to start with as I need to change a load of other code:
  return( al->left_scores[row] + gap_penalty( al->left_tr_gap[row] + 1, (float)pen->min_intron_size, 
					      pen->gap_i, pen->gap_e_mod, pen->fshift) );
}

float score_from_above(struct al_info *al, struct al_penalty *pen, int row, int g_i){
  return( al->right_scores[row-1] + 
	  gap_penalty( al->right_gen_gap[row-1] + 1, (float)pen->min_intron_size, 
		       pen->gap_i, pen->gap_e_mod, pen->fshift) +
	  horizontal_gap_close_penalty(al, pen, row, al->right_tr_gap, g_i)
	  );
}

// This will need to do the same checks for gap sequences as performed in score_from_left
// and score_from_above. This suggest that I should refactor that code into their own
// functions. Eg. check_vert_gap, check_horizontal_gap
float score_from_above_left(struct al_info *al, struct al_penalty *pen, int row, int g_i){
  // First get the substitution matrix score. Note that the we don't have to check for
  // this row being a gap one as the substitution matrix will 
  unsigned int tr_i = row - 1;
  int m_score_o = (al->tr[tr_i] - pen->al_offset) + (al->gen[g_i] - pen->al_offset) * pen->al_size;
  float diag_score = al->left_scores[row-1] + (float)pen->sub_table[ m_score_o ];
  // add on penalties for vertical or gap_close_penalties; we could check if we need to do this before
  // we call the functions to avoid the function call, but first let's just do it:
  return( diag_score + // vertical_gap_close_penalty( al, pen, row-1 ) +
	  horizontal_gap_close_penalty( al, pen, row, al->left_tr_gap, g_i-1 ));
}
 
// Like align_transcript_bitp, but splice site aware
// Considers intron insertion to be free only when GT-AG, AT-AC, GC-AG
// 
// This uses gap-closing score calculation; that is gap penalties are
// only calculated upon the closing of the gap. This will make this
// implementation heuristic. To guarantee optimal scores we would need
// to look back one more step, and that would still be problematic.
// if global is non-0, then do a semi-global alignment (row 0 is 0,
// col 0 is set according to gap extension rules, and alignment is taken from
// the highest score in the last row).
struct sw_alignment *align_transcript_bitp_ssa(const unsigned char *tr, const unsigned char *gen,
					       int tr_l, int gen_l, float gap_i, float gap_e_mod, float fshift,
					       float intron_new, float intron_bad, int min_intron_size,
					       int *sub_table, int al_offset, int al_size, char intron_char,
					       int global)
{
  struct al_info al = init_al_info(tr, gen, tr_l, gen_l, intron_char);
  const unsigned char intron_type_n = INTRON_TYPE_N;
  const char *intron_start[INTRON_TYPE_N] = {"GT", "AT", "GC"};
  const char *intron_end[INTRON_TYPE_N] = {"AG", "AC", "AG"};
  struct al_penalty pen = init_al_penalty(sub_table, al_offset, al_size, intron_char,
					  gap_i, gap_e_mod, fshift, intron_new, intron_bad, min_intron_size,
					  intron_type_n, intron_start, intron_end
					  );
  // the pointers.. 
  struct bit_ptr ptr = init_bit_ptr( tr_l + 1, gen_l + 1);
  // if global prepare the first column:
  if(global){
    for(int i=1; i <= tr_l; ++i){
      al.left_scores[i] = score_from_above(&al, &pen, i, 0);
      // note that score_from_above will call horizontal_gap_close_penalty with
      // right_tr_gap as an argument. This is WRONG, but doesn't matter. However,
      // future changes could cause a problem here.
      al.left_gen_gap[i] = al.left_gen_gap[i-1] + 1;
      set_ptr(&ptr, i, 0, (unsigned char)2);
    }
  }
  // and the maximum position;
  float max_score = 0;
  int max_col = 0;
  int max_row = 0;
  // if doing a (semi) global analysis we need to find the max in the last row:
  float last_row_max = al.left_scores[tr_l];
  int last_row_max_col = 0;
  // if we are doing a local alignment, then we need to include a 0 alternative:
  int first_score_alt = global ? 1 : 0;
  // This should be enough...
  for(int g_i=0; g_i < gen_l; ++g_i){
    int col = g_i + 1;
    for(int tr_i=0; tr_i < tr_l; ++tr_i){
      int row = tr_i + 1;
      // left_count_prev[row] = left_count[row];
      // m_score_o, is the offset in the substitution matrix
      // int m_score_o = (tr[tr_i] - al_offset) + (gen[g_i] - al_offset) * al_size;
      
      float left_score = score_from_left(&al, &pen, row);
      float up_score = score_from_above(&al, &pen, row, g_i);
      float diag_score = score_from_above_left(&al, &pen, row, g_i);

      float scores[4] = {0, left_score, up_score, diag_score};
      int max_i = first_score_alt;
      for(int i=first_score_alt + 1; i < 4; ++i)
	max_i = (scores[i] > scores[max_i]) ? i : max_i;
      // If 3 or 0, then both gap counts will be set to 0.
      al.right_tr_gap[row] = (max_i == 1) ? al.left_tr_gap[row] + 1 : 0;
      al.right_gen_gap[row] = (max_i == 2) ? al.right_gen_gap[row-1] + 1 : 0;
      // Check if we have crossed any intron characters;
      if(max_i == 2){
	if(tr[tr_i] == intron_char)
	  al.right_int_x[row]++;
      }else{
	al.right_int_x[row] = 0;
      }
      // and set the max score and the pointer
      al.right_scores[row] = scores[max_i];
      set_ptr( &ptr, row, col, (unsigned char)max_i );
      // update the maximum if appropriate;
      if(al.right_scores[row] > max_score){
	max_score = al.right_scores[row];
	max_col = col;
	max_row = row;
      }
    }
    if(al.right_scores[tr_l] > last_row_max){
      last_row_max = al.right_scores[tr_l];
      last_row_max_col = col;
    }
    swap_left_right( &al );
  }
  
  struct sw_alignment *alignment = malloc(sizeof(struct sw_alignment));
  memset(alignment, 0, sizeof(struct sw_alignment));
  int error = 0;
  if(!global)
    error = extract_bitp_alignment( tr, gen, &ptr, max_row, max_col, max_score, alignment);
  else
    error = extract_bitp_alignment( tr, gen, &ptr, tr_l, last_row_max_col, max_score, alignment);
  // clear all the allocated memory:
  free_al_info(&al);
  free_al_penalty(&pen);
  free_bit_ptr(&ptr);
  return(alignment);
}
