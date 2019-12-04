#include <R.h>
#include <Rinternals.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include "needleman_wunsch.h"

// This is a rewrite of exon_aligneR.c
// Modified to:
// 1. Use filly by column matrices like in R
// 2. Always place seq 1 (a) by rows and seq 2 (b) by columns
// 3. Remove the merge functionaliy


struct transcript {
  const char *seq;
  int t_length;
  int e_n;
  int *e_lengths;
  int *e_offsets;
};

struct transcript transcript_init(SEXP seq, SEXP e_lengths){
  if(TYPEOF(seq) != STRSXP || length(seq) != 1)
    error("seq should be of type STRSXP and have a length of 1");
  if(TYPEOF(e_lengths) != INTSXP || length(e_lengths) < 1)
    error("e_lengths must be an integer vector of length > 0");
  struct transcript trans;
  SEXP s = STRING_ELT(seq, 0);
  trans.t_length = length(s);
  trans.seq = CHAR(s);
  trans.e_n = length( e_lengths );
  trans.e_lengths = INTEGER( e_lengths );
  // check that sum is equal to t_length;
  int e_sum = 0;
  for(int i=0; i < trans.e_n; ++i)
    e_sum += trans.e_lengths[i];
  if(e_sum != trans.t_length)
    error("exon lengths do no add up to total transcript length");
  trans.e_offsets = malloc( sizeof(int) * trans.e_n );
  trans.e_offsets[0] = 0;
  for(int i=1; i < trans.e_n; ++i)
    trans.e_offsets[i] = trans.e_offsets[i-1] + trans.e_lengths[i-1];
  return( trans );
};

void transcript_free(struct transcript trans){
  free(trans.e_offsets);
}

char* mk_exon(struct transcript trans, int i){
  int l = trans.e_lengths[i];
  char *exon = malloc(sizeof(char) * (1 + l));
  exon[l] = 0;
  memcpy((void*) exon, (void*)(trans.seq + trans.e_offsets[i]), sizeof(char) * l );
  return(exon);
}

char* mk_gaps(int g_n){
  char *exon = malloc(sizeof(char) * (g_n + 1));
  exon[g_n] = 0;
  memset((void*)exon, '-', g_n);
  return(exon);
}

// al_length contains the length of the alignment
// *alignment is a matrix containing the exon indices
// at the different positions of the alignment. A gap
// is indicated by -1.
// the alignment is filled by column as in R.
struct gene_alignment {
  int length;
  int *alignment;
};

struct gene_alignment gene_alignment_init(){
  struct gene_alignment ga;
  ga.length = 0;
  ga.alignment = 0;
  return(ga);
}

struct aligned_exons {
  int length;
  int *lengths;
  char **a;
  char **b;
};

struct aligned_exons aligned_exons_init(int n){
  struct aligned_exons al;
  al.length = n;
  al.lengths = malloc(sizeof(int) * n);
  al.a = malloc(sizeof(char*) * n);
  al.b = malloc(sizeof(char*) * n);
  return(al);
}

void aligned_exons_free(struct aligned_exons al){
  for(int i=0; i < al.length; ++i){
    free(al.a[i]);
    free(al.b[i]);
  }
  free(al.lengths);
  free(al.a);
  free(al.b);
}


// takes pre-allocated tables
// take separate values for horizontal and vertical gap penalties to
// allow funny things.. 
void init_nm_tables( double *scores, char *pointers, int m_height, int m_width,
		     double h_gap_i, double h_gap_e, double v_gap_i, double v_gap_e,
		     const char left, const char up){
  memset( (void*)scores, 0, sizeof(double) * m_height * m_width );
  memset( (void*)pointers, 0, sizeof(char) * m_height * m_width);

  // Set the initial scores and pointers
  scores[ m_offset(0, 1, m_height) ] = h_gap_i;
  scores[ m_offset(1, 0, m_height) ] = v_gap_i;
  pointers[ m_offset(0, 1, m_height) ] = left;
  pointers[ m_offset(1, 0, m_height) ] = up;

  // and then fill in the first row and column
  for(int i=2; i < m_width; ++i){
    scores[ m_offset(0, i, m_height) ] = scores[ m_offset(0, i-1, m_height) ] + h_gap_e;
    pointers[ m_offset(0, i, m_height) ] = left;
  }
  for(int i=2; i < m_height; ++i){
    scores[ m_offset(i, 0, m_height) ] = scores[ m_offset(i-1, 0, m_height) ] + v_gap_e;
    pointers[ m_offset(i, 0, m_height) ] = up;
  }

}

// dp for dynamic programming
struct dp_max {
  int row;
  int column;
  double max;
};

struct dp_max dp_max_init(){
  struct dp_max max;
  max.row=0;
  max.column=0;
  max.max=0;
  return(max);
}

void dp_max_update( struct dp_max *max, double score, int row, int column ){
  if(max->max <= score){
    max->max = score;
    max->row = row;
    max->column = column;
  }
}

// expects to be given a table for the pointers since this is necessary to trace the alignment back
// if local is 1, the alignment is a localised nm where terminal gaps are not penalised. I don't
// know the proper name of such an alignment, but it may be suitable here..
struct dp_max exon_nm(const char *seq_a, const char *seq_b, int a_l, int b_l, double *penalties,
	       double *scores, char *pointers, char local){
  // the penalty values
  double match = penalties[0];
  double mis_match = penalties[1];
  double gap_i = penalties[2];
  double gap_e = penalties[3];
  const char left = 1;
  const char up = 2;
  // the dimensions of the table
  int m_height = a_l + 1;
  int m_width = b_l + 1;
  
  // use a full score table. Do not try to minimise memory here. Keep things simply.
  // We assume that both the pointers and the scores tables have been set up with proper dimensions.
  if(!local){
    init_nm_tables( scores, pointers, m_height, m_width, gap_i, gap_e, gap_i, gap_e, left, up);
  }else{
    double h_gap_i=gap_i, h_gap_e=gap_e, v_gap_i=gap_i, v_gap_e=gap_e;
    if(m_width > m_height){
      h_gap_i = 0;
      h_gap_e = 0;
    }
    if(m_height > m_width){
      v_gap_i = 0;
      v_gap_e = 0;
    }
    init_nm_tables( scores, pointers, m_height, m_width, h_gap_i, h_gap_e, v_gap_i, v_gap_e, left, up);
  }
  // And then we simply go through the table positions. Let us make a pointer to the two
  // sequences that we are using
  int o=0, o_l, o_u, o_d;  // offsets for the different positions
  struct dp_max max = dp_max_init();
  double left_gap=0, up_gap=0;
  
  for(int row=1; row < m_height; ++row){
    for(int column=1; column < m_width; ++column){
      double sc[3];
      o = m_offset(row, column, m_height);
      o_l = m_offset(row, column-1, m_height);
      o_u = m_offset(row-1, column, m_height);
      o_d = m_offset(row-1, column-1, m_height);
      // We only want to allow gaps for one of the directions. I.e. the minimal score should not
      // be 0. 
      left_gap = (m_width > m_height && local && row == m_height - 1) ? 0 : ( pointers[o_l] == left ? gap_e : gap_i );
      up_gap = (m_height > m_width && local && column == m_width - 1) ? 0 : ( pointers[o_u] == up ? gap_e : gap_i );
      
      sc[0] = scores[o_l] + left_gap; // ( pointers[o_l] == left ? gap_e : gap_i );
      sc[1] = scores[o_u] + up_gap; // ( pointers[o_u] == up ? gap_e : gap_i );
      sc[2] = scores[o_d] + (seq_a[row-1] == seq_b[column-1] ? match : mis_match );
      int max_i = which_max_d(sc, 3);
      scores[o] = sc[max_i];
      pointers[o] = max_i + 1;
    }
    //    dp_max_update( &max, scores[o], row, m_width - 1 );
  }
  //  if(!local){
  max.max = scores[o];
  max.row = m_height - 1;
  max.column = m_width - 1;
  return(max);
}

// the use of int for the pointers is wasteful, but it makes it easier to
// interface with R. I would otherwise have to build a table of char*,
// which is even worse..
// if local == TRUE then perform a Smith Waterman
struct dp_max gene_align( struct transcript a, struct transcript b, double *exon_scores,
		 double match, double gap,
		 double *scores, int *pointers, char local ){
  int height = a.e_n + 1;
  int width = b.e_n + 1;
  int m_size = width * height;
  memset( (void*)scores, 0, sizeof(double) * m_size );
  memset( (void*)pointers, 0, sizeof(int) * m_size );

  int left = 1;
  int up = 2;
  // init the pointers.. if not local
  if(!local){
    for(int row=1; row < height; ++row){
      int o = m_offset( row, 0, height );
      int o_u = m_offset( row-1, 0, height );
      pointers[o] = up;
      scores[o] = scores[o_u] - match * a.e_lengths[row-1] * gap;
    }
    for(int column=1; column < width; ++column){
      int o = m_offset( 0, column, height );
      int o_l = m_offset( 0, column - 1, height );
      pointers[o] = left;
      scores[o] = scores[o_l] - match * b.e_lengths[column-1] * gap;
    }
  }
  int o=0;
  struct dp_max max_pos = dp_max_init();
  for(int row=1; row < height; ++row){
    for(int column=1; column < width; ++column){
      o = m_offset(row, column, height);
      int o_l = m_offset(row, column-1, height);
      int o_u = m_offset(row-1, column, height);
      int o_d = m_offset(row-1, column-1, height);
      double sc[4];
      sc[0] = 0;
      sc[1] = scores[o_l] - match * b.e_lengths[column-1] * gap;
      sc[2] = scores[o_u] - match * a.e_lengths[row-1] * gap;
      sc[3] = scores[o_d] + exon_scores[ m_offset(row-1, column-1, height-1) ];
      /* Rprintf("row: %d column %d previous scores: %f, %f, %f\n", row, column, scores[o_l], scores[o_u], scores[o_d]); */
      /* Rprintf("sc: %f, %f, %f, %f\n", sc[0], sc[1], sc[2], sc[3]); */
      if(!local){
	int max_i = which_max_d( &sc[1], 3 );
	scores[o] = sc[ max_i +1 ];
	pointers[o] = max_i + 1;
	//	Rprintf("not local, max_i %d  max: %f\n", max_i, sc[ max_i + 1 ]);
      }else{
	int max_i = which_max_d( sc, 4 );
	scores[o] = sc[ max_i ];
	pointers[o] = max_i;
	//	Rprintf("local, max_i %d  max: %f\n", max_i, sc[ max_i ]);
      }
      dp_max_update( &max_pos, scores[o], row, column );
    }
  }
  // override the global position if not global
  if(!local){
    max_pos.max = scores[o];
    max_pos.row = height - 1;
    max_pos.column = width - 1;
  }
  return( max_pos );
}



// height and width are the dimensions of the table, not the number of exons
void extract_gene_alignment(int *pointers, int height, int width, struct gene_alignment* align, struct dp_max max_pos){
  int al_length = 0;
  /* int row = height -1; */
  /* int column = width -1; */
  int row = max_pos.row;
  int column = max_pos.column;
  while(row > 0 || column > 0){
    int o = m_offset( row, column, height );
    if(pointers[o] == 0)
      break;
    row = (pointers[o] &  2) ? row - 1 : row;
    column = (pointers[o] & 1) ? column - 1 : column;
    al_length++;
  }
  if(!al_length)
    return;
  align->length = al_length;
  align->alignment = malloc( sizeof(int) * 2 * al_length );
  int *al_a = align->alignment;
  int *al_b = align->alignment + al_length;
  row = max_pos.row;
  column = max_pos.column;
  /* row = height -1; */
  /* column = width -1; */
  while(row > 0 || column > 0){
    al_length--;
    int o = m_offset( row, column, height );
    if(pointers[o] == 0 || al_length < 0)
      break;
    // This asks the same question twice and is hence bad. But the
    // only alternative I can think of makes use of two if-else
    // constructs and that too is ugly.
    al_a[al_length] = (pointers[o] & 2) ? row - 1 : -1;
    al_b[al_length] = (pointers[o] & 1) ? column - 1 : -1;
    row = (pointers[o] &  2) ? row - 1 : row;
    column = (pointers[o] & 1) ? column - 1 : column;
  }
}

// It should be possible to use for both local and global alignments.
void extract_exon_alignment(char *pointers, int height, int width, const char *a, const char *b, struct dp_max max,
			    char **a_a, char **b_a){
  int al_length = 0;
  int row = max.row;
  int column = max.column;
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
  row = max.row;
  column = max.column;
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

// Assumes that we never get -1 and -1 for both exons. That would be an error that we should
// check for. But let's write some working code first and then harden it.
struct aligned_exons extract_exon_alignments(struct transcript a, struct transcript b, struct gene_alignment g_align, double *penalties, char local){
  // For each row of the gene_alignment table we want to create two char* structs
  // containing the aligned sequences.
  struct aligned_exons aligns = aligned_exons_init( g_align.length );
  int *ex_a = g_align.alignment;
  int *ex_b = g_align.alignment + g_align.length;
  double *scores = 0;
  char *pointers = 0;
  for(int i=0; i < g_align.length; ++i){
    int a_i = ex_a[i];
    int b_i = ex_b[i];
    if(a_i == -1){
      aligns.lengths[i] = b.e_lengths[b_i];
      aligns.a[i] = mk_gaps( b.e_lengths[ b_i ] );
      aligns.b[i] = mk_exon( b, b_i );
      continue;
    }
    if(b_i == -1){
      aligns.lengths[i] = a.e_lengths[a_i];
      aligns.a[i] = mk_exon( a, a_i );
      aligns.b[i] = mk_gaps( a.e_lengths[a_i] );
      continue;
    }
    // Here I need to align the two exons again.
    int a_l = a.e_lengths[ a_i ];
    int b_l = b.e_lengths[ b_i ];
    int m = (a_l + 1) * (b_l + 1);
    scores = realloc( scores, sizeof(double) * m );
    pointers = realloc( pointers, sizeof(char) * m);
    struct dp_max max = exon_nm( (a.seq + a.e_offsets[a_i]), (b.seq + b.e_offsets[b_i]),
	     a.e_lengths[a_i], b.e_lengths[b_i], penalties, scores, pointers, local );
    // then traverse the score table and work out the length of the alignment...
    extract_exon_alignment(pointers, a_l + 1, b_l + 1, a.seq + a.e_offsets[a_i], b.seq + b.e_offsets[b_i], max,
			   aligns.a + i, aligns.b + i);
  }
  free(scores);
  free(pointers);
  return(aligns);
}

// a_seq_r and b_seq_r: transcript sequences
// a_lengths, b_lengths : the length of exons in transcripts a and b
// e_penalties_r: the penalties / scores used for aligning sequences to each other
// g_penalties_r: some parameters that can be used to define the penalties
//                for the final exon alignments
SEXP align_exons(SEXP a_seq_r, SEXP b_seq_r, SEXP a_lengths, SEXP b_lengths,
		 SEXP e_penalties_r, SEXP g_penalty_r, SEXP local_r, SEXP local_gene_r ){
  // SANITY check
  if( TYPEOF(e_penalties_r) != REALSXP || TYPEOF(g_penalty_r) != REALSXP )
    error("Arguments 3 and 4 should vectors of real numbers");
  if(length(e_penalties_r) != 4)
    error("The third argument should provide: match, mismatch, gap insertion, gap extension values");
  if(length(g_penalty_r) != 1)
    error("The fourth argument should provide a single value from which we can derive an exon gap penalty");
  if(TYPEOF(local_r) != LGLSXP)
    error("The fifth argument should be a logical value determining whether we use localised NM-alignments");
  if(TYPEOF(local_gene_r) != LGLSXP)
    error("The sixth argument should be a logical value determining whether we use localised NM-alignments");
  
  // the transcript_init function sanity checks its SEXP arguments
  struct transcript trans_a = transcript_init(a_seq_r, a_lengths);
  struct transcript trans_b = transcript_init(b_seq_r, b_lengths);
  double *e_penalties = REAL(e_penalties_r);
  double g_penalty = REAL(g_penalty_r)[0];
  char local = (char)asLogical(local_r);
  char local_gene = (char)asLogical(local_gene_r);
  
  // Allocate space for a table of scores that we can return to R
  int a_n = trans_a.e_n;
  int b_n = trans_b.e_n;
  SEXP ret_data = PROTECT( allocVector( VECSXP, 6 ));
  SET_VECTOR_ELT(ret_data, 0, allocMatrix(REALSXP, a_n, b_n ));
  SET_VECTOR_ELT(ret_data, 1, allocMatrix(REALSXP, a_n + 1, b_n + 1 ));
  SET_VECTOR_ELT(ret_data, 2, allocMatrix(INTSXP, a_n + 1, b_n + 1 ));

  // Get the resulting matrix:
  double *exon_scores = REAL( VECTOR_ELT(ret_data, 0));
  double *gene_score_matrix = REAL( VECTOR_ELT(ret_data, 1));
  int *gene_pointer_matrix = INTEGER( VECTOR_ELT(ret_data, 2));
  double *scores = 0;
  char *pointers = 0;
  
  for(int row=0; row < a_n; ++row){
    for(int column=0; column < b_n; ++column){
      //      int m_size = (1 + a_exons.seq_l[row]) * (1 + b_exons.seq_l[column]);
      int m_size = (1 + trans_a.e_lengths[row]) * (1 + trans_b.e_lengths[column]);
      scores = realloc((void*)scores, sizeof(double) * m_size );
      pointers = realloc((void*)pointers, sizeof(char) * m_size );
      struct dp_max max  = exon_nm( trans_a.seq + trans_a.e_offsets[row],
				    trans_b.seq + trans_b.e_offsets[column],
				    trans_a.e_lengths[row], trans_b.e_lengths[column],
				    e_penalties, scores, pointers, local );
      exon_scores[ m_offset(row, column, a_n) ] = max.max;
    }
  }

  struct dp_max max_pos = gene_align( trans_a, trans_b, exon_scores, e_penalties[0], g_penalty, gene_score_matrix, gene_pointer_matrix, local_gene );
  struct gene_alignment g_align = gene_alignment_init();
  extract_gene_alignment( gene_pointer_matrix, trans_a.e_n + 1, trans_b.e_n + 1, &g_align, max_pos );
  SET_VECTOR_ELT(ret_data, 3, allocMatrix( INTSXP, g_align.length, 2 ));
  memcpy( (void*)INTEGER(VECTOR_ELT(ret_data, 3)), (void*)g_align.alignment, sizeof(int) * 2 * g_align.length );

  struct aligned_exons al_exons = extract_exon_alignments(trans_a, trans_b, g_align, e_penalties, local);
  SET_VECTOR_ELT(ret_data, 4, allocVector(STRSXP, al_exons.length) );
  SET_VECTOR_ELT(ret_data, 5, allocVector(STRSXP, al_exons.length) );
  for(int i=0; i < al_exons.length; ++i){
    SET_STRING_ELT( VECTOR_ELT(ret_data, 4), i, mkChar( al_exons.a[i] ));
    SET_STRING_ELT( VECTOR_ELT(ret_data, 5), i, mkChar( al_exons.b[i] ));
  }
  
  aligned_exons_free(al_exons);
  
  free(scores);
  free(pointers);
  free(g_align.alignment);
  transcript_free( trans_a );
  transcript_free( trans_b );
  UNPROTECT(1);
  return(ret_data);
}

SEXP align_seqs(SEXP a_seq_r, SEXP b_seq_r, SEXP al_offset_r, SEXP al_size_r,
		SEXP sub_matrix_r, SEXP gap_r,
		SEXP tgaps_free_r, SEXP special_char_r){
  if( TYPEOF(a_seq_r) != STRSXP || TYPEOF(b_seq_r) != STRSXP || TYPEOF(special_char_r) != STRSXP )
    error("a_seq_r, b_seq_r and special_char_r must all be R strings");
  if( TYPEOF(al_offset_r) != INTSXP || TYPEOF(al_size_r) != INTSXP )
    error("al_offset_r and al_size_r must both be ints");
  if( TYPEOF(sub_matrix_r) != INTSXP || !isMatrix(sub_matrix_r) )
    error("sub_matrix_r should be a numeric matrix");
  if( TYPEOF(gap_r) != INTSXP || length(gap_r) != 2 )
    error("gap_r should be an integer vector of length 2 (insertion, extension)");
  if( TYPEOF(tgaps_free_r) != LGLSXP )
    error("tgaps_free_r should be a logical vector");
  if(length(al_offset_r) != 1 || length(al_size_r) != 1 || length(a_seq_r) != 1
     || length(b_seq_r) != 1 || length(tgaps_free_r) != 1 || length(special_char_r) != 1 )
    error("most arguments should have a length of one");

  int al_offset = asInteger(al_offset_r);
  int al_size = asInteger(al_size_r);
  
  SEXP sub_matrix_dims = getAttrib(sub_matrix_r, R_DimSymbol);
  int nrow = INTEGER(sub_matrix_dims)[0];
  int ncol = INTEGER(sub_matrix_dims)[1];
  if(nrow != ncol || nrow != al_size)
    error("faulty matrix dimensions: %d, %d should both be: %d", nrow, ncol, al_size);
  int *sub_table = INTEGER( sub_matrix_r );
  int tgaps_free = asLogical( tgaps_free_r );

  int a_l = length( STRING_ELT( a_seq_r, 0 ));
  const unsigned char *a_seq = (const unsigned char*)CHAR( STRING_ELT( a_seq_r, 0 ));
  int b_l = length( STRING_ELT( b_seq_r, 0 ));
  const unsigned char *b_seq = (const unsigned char*)CHAR( STRING_ELT( b_seq_r, 0 ));

  int *gap = INTEGER(gap_r);
  
  // We can then allocate the required memorey. No chars everything is in ints
  SEXP ret_data = PROTECT(allocVector( VECSXP, 7 ));
  SET_VECTOR_ELT( ret_data, 0, allocMatrix( INTSXP, a_l+1, b_l+1 ));
  SET_VECTOR_ELT( ret_data, 1, allocMatrix( INTSXP, a_l+1, b_l+1 ));
  SET_VECTOR_ELT( ret_data, 2, allocVector( STRSXP, 2 ));

  int *score_table = INTEGER( VECTOR_ELT( ret_data, 0 ));
  int *ptr_table = INTEGER( VECTOR_ELT( ret_data, 1 ));

  needleman_wunsch( a_seq, b_seq, a_l, b_l, gap[0], gap[1], sub_table, al_offset, al_size, tgaps_free,
  		    score_table, ptr_table );

  char *al, *bl;  // a and b aligned
  struct align_stats al_stats = extract_nm_alignment( ptr_table, a_l+1, b_l+1, (const char*)a_seq, (const char*)b_seq, &al, &bl );

  SET_STRING_ELT( VECTOR_ELT(ret_data, 2), 0, mkChar( al ));
  SET_STRING_ELT( VECTOR_ELT(ret_data, 2), 1, mkChar( bl ));

  int special_n = length( STRING_ELT( special_char_r, 0 ) );
  if( special_n > 0 ){
    const char *special_char = CHAR( STRING_ELT( special_char_r, 0 ) );
    int a_sl, b_sl;
    int *a_s, *b_s;
    char_at( al, special_char[0], &a_s, &a_sl );
    char_at( bl, special_char[0], &b_s, &b_sl );
    int al_nrows;
    int *al_table = aligned_i( a_s, b_s, a_sl, b_sl, &al_nrows );
    SET_VECTOR_ELT( ret_data, 3, allocVector( INTSXP, a_sl ));
    SET_VECTOR_ELT( ret_data, 4, allocVector( INTSXP, b_sl ));
    SET_VECTOR_ELT( ret_data, 5, allocMatrix( INTSXP, al_nrows, 2 ));
    memcpy( (void*)INTEGER(VECTOR_ELT(ret_data, 3)), (void*)a_s, sizeof(int) * a_sl);
    memcpy( (void*)INTEGER(VECTOR_ELT(ret_data, 4)), (void*)b_s, sizeof(int) * b_sl);
    memcpy( (void*)INTEGER(VECTOR_ELT(ret_data, 5)), (void*)al_table, sizeof(int) * 2 * al_nrows );
    free(a_s);
    free(b_s);
    free( al_table );
  }

  // and then let us get the stats
  SET_VECTOR_ELT( ret_data, 6, allocVector( INTSXP, sizeof(al_stats) / sizeof(int) ));
  memcpy( (void*)INTEGER(VECTOR_ELT(ret_data, 6)), &al_stats, sizeof(al_stats) );

  free(al);
  free(bl);
  
  UNPROTECT(1);
  return( ret_data);
}

// To make a mulithreaded version I need to first define a struct that takes the arguments
// that need to be passed to needlmean_wunsch and extract_nm_alignment
// Note that since this is for high-throughput alignments we do not return the scoring matrices
// so we do not need to have these passed in as arguments

// we do not really need an initialiser as we can simply use memset here
struct needle_wunsch_results {
  char *al_a, *al_b;      // sequences with gaps
  int *a_pos, *b_pos;     // positions of special characters
  int a_pos_l, b_pos_l;   // number of special characters in each sequence
  int *pos_table;         // a table of the ones that have been aligned.
  int pos_table_nrow;     // the number of rows in the table
  int score;
  struct align_stats al_stats;
};

struct needle_wunsch_args {
  // to create the scoring matrix:
  const unsigned char *a, **b;    // the sequences
  int a_l;                       // length of the single sequence a
  int b_n;                       // number of b sequences
  int *b_l;                      // the lengths of the b sequences
  int gap_i, gap_e;              // gap insertion and gap extension penalties
  int *sub_table;                // the substitution matrix
  int al_offset, al_size;        // the alphabet offset and size
  int tgaps_free;                // 0/1: should terminal gaps in the shorter sequence be free
  char special_char;
  // to extract the alignment:
  struct *needle_wunsch_results ret_data;  // b_n results
  unsigned char *jobs;                      // b_n jobs, 0 if not yet done
  pthread_mutex_t *mutex;
};



// 
struct needle_wunsch_args needle_wunsch_args_init(const unsigned char *a, int a_l, const unsigned char **b, int b_n, int *b_l,
						    int gap_i, int gap_e, int *sub_table, int al_offset, int al_size, int tgaps_free,
						    char special_char, pthread_mutex_t *mutex
						    ){
  struct needle_wunsch_args args;
  args.a = a;
  args.a_l = a_l;
  args.b = b;
  args.b_n = b_n;
  args.b_l = b_l;
  args.gap_i = gap_i;
  args.gap_e = gap_e;
  args.sub_table = sub_table;
  args.al_offset = al_offset;
  args.al_size = al_size;
  args.tgaps_free = tgaps_free;
  args.special_char = special_char;
  args.mutex = mutex;
  
  // other pointers are best initialised to 0, so that freeing them doesn't hurt..
  args.ret_data = malloc( sizeof(struct needle_wunsch_results) * b_n );
  args.jobs = malloc( sizeof(unsigned char) * b_n );
  // and initialise the data to something useful.. 
  memset( (void*)args.ret_data, 0, sizeof(struct needle_wunsch_results) * b_n );
  memset( (void*)args.jobs, 0, sizeof(unsigned char) * b_n );

}

void needle_wunsch_args_free( struct *needle_wunsch_args args ){
  for(int i=0; i < args.b_n; ++i){
    free( args->ret_data[i].al_a )
    free( args->ret_data[i].al_b )
    free( args->ret_data[i].a_pos )
    free( args->ret_data[i].b_pos )
    free( args->ret_data[i].pos_table )
  }
  free( args->ret_data );
  free( args->jobs );
}

// And a method that performs a single alignment for the given sequences..
void * needleman_wunsch_thread(void *args_ptr){
  struct needle_wunsh_args *args = (struct needle_wunsh_args*)args_ptr;
  int score_table = 0;
  int ptr_table = 0;
  int i = 0; // the current job
  while(1){
    pthread_mutex_lock( args->mutex );
    while( i < args->b_n && args->jobs[ i ] )
      ++i;
    if(args->jobs[i])
      ++i;             // we can check..
    else
      args->jobs[i] = 1;
    pthread_mutex_unlock( args->mutex );
    if( i >= args->b_n )
      break;
    
    int height = args->a_l + 1;
    int width = args->b_l[i] + 1;
    *score_table = realloc((void*)scrore_table, sizeof(int) * width * height );
    *ptr_table = realloc((void*)ptr_table, sizeof(int) * width * height );
    
    needleman_wunsch( args->a, args->b[i], args->a_l, args->b_l[i],
		      args->gap_i, args->gap_e,
		      args->sub_table, args->al_offset, args->al_size,
		      score_table, ptr_table);
    
    struct needle_wunsch_results *ret_data = &(args->ret_data[i]);
    ret_data->score = score_table[ width * height - 1 ];
    ret_data->al_stats = extract_nm_alignment(ptr_table, height, width, args->a, args->b[i],
					      &(ret_data->al_a), &(ret_data->al_b));
    
    char_at( ret_data->al_a, args->special_char, &(ret_data->pos_a), &(ret_data->pos_a_l));
    char_at( ret_data->al_b, args->special_char, &(ret_data->pos_b), &(ret_data->pos_b_l));
    ret_data->pos_table = aligned_i( ret_data->pos_a, ret_data->pos_b, ret_data->pos_a_l, ret_data->pos_b_l,
				     &(ret_data->pos_table_nrow) );
  }
  free(score_table);
  free(ptr_table);
  pthread_cond_signal( &cond );
  pthread_exit((void*)args_ptr);
}
    

SEXP align_seqs_mt(SEXP a_seq_r, SEXP b_seq_r, SEXP al_offset_r, SEXP al_size_r,
		   SEXP sub_matrix_r, SEXP gap_r,
		   SEXP tgaps_free_r, SEXP special_char_r, SEXP n_thread_r){
  if( TYPEOF(a_seq_r) != STRSXP || TYPEOF(b_seq_r) != STRSXP || TYPEOF(special_char_r) != STRSXP )
    error("a_seq_r, b_seq_r and special_char_r must all be R strings");
  if( TYPEOF(al_offset_r) != INTSXP || TYPEOF(al_size_r) != INTSXP )
    error("al_offset_r and al_size_r must both be ints");
  if( TYPEOF(sub_matrix_r) != INTSXP || !isMatrix(sub_matrix_r) )
    error("sub_matrix_r should be a numeric matrix");
  if( TYPEOF(gap_r) != INTSXP || length(gap_r) != 2 )
    error("gap_r should be an integer vector of length 2 (insertion, extension)");
  if( TYPEOF(tgaps_free_r) != LGLSXP )
    error("tgaps_free_r should be a logical vector");
  if( TYPEOF(n_thread_r) != INTSXP || length(n_thread_r) != 1)
    error("n_thread_r should be an integer vector of length 1");
  if(length(al_offset_r) != 1 || length(al_size_r) != 1 || length(a_seq_r) != 1
     || length(b_seq_r) < 1 || length(tgaps_free_r) != 1 || length(special_char_r) != 1 )
    error("most arguments should have a length of one (b_seq_r may have more than one sequence)");

  int al_offset = asInteger(al_offset_r);
  int al_size = asInteger(al_size_r);
  
  SEXP sub_matrix_dims = getAttrib(sub_matrix_r, R_DimSymbol);
  int nrow = INTEGER(sub_matrix_dims)[0];
  int ncol = INTEGER(sub_matrix_dims)[1];
  if(nrow != ncol || nrow != al_size)
    error("faulty matrix dimensions: %d, %d should both be: %d", nrow, ncol, al_size);
  int *sub_table = INTEGER( sub_matrix_r );
  int tgaps_free = asLogical( tgaps_free_r );

  if(!length(STRING_ELT(special_char_r, 0)))
    error("A special character mush be given");

  char special_char = CHAR( STRING_ELT( special_char_r, 0) )[0];
  int a_l = length( STRING_ELT( a_seq_r, 0 ));
  const unsigned char *a_seq = (const unsigned char*)CHAR( STRING_ELT( a_seq_r, 0 ));

  int *gap = INTEGER(gap_r);
  int n_thread = asInteger( n_thread_r );
  if(n_thread <= 1)
    error("n_thread must be larger than 0!");
  
  // We now need to set up vectors containing the appropriate b sequences
  int b_n = length(b_seq_r);
  int *b_l = malloc( sizeof(int) * b_n );
  const unsigned char **b_seq = malloc( sizeof(const char*) * b_n );
  for(int i=0; i < b_n; ++i){
    b_l[i] = length( STRING_ELT( b_seq_r, i ));
    b_seq[i] = (const unsigned char*)CHAR( STRING_ELT( b_seq_r, 0 ));
  }
  // reduce the thread number if we have fewer b_sequences than threads
  if( b_n < n_thread )
    n_thread = b_n;

  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  // The we have to set up arrays of pthreads and needle_wunsch_args
  struct needle_wunsch_args args = needle_wunsch_args_init(a_seq, a_l, b_seq, b_n, b_l, gap[0], gap[1],
							   sub_table, al_offset, al_size, tgaps_free,
							   special_char, &mutex);
							   
  pthread_t *threads = malloc( sizeof(pthread_t), n_thread );
  for(int i=0; i < n_thread; ++i){
    pthread_create( &threads[i], NULL, needleman_wunsch_thread, (void*)&args );
  }
  // then we wait for all threads t exit..
  void *status;
  for(int i=0; i < n_thread; ++i)
    pthread_join(threads[i], &status);

  // 1. Create the R data structures
  // 2. Copy over the data to these
  // 3. Clean up allocated data and then return.
  SEXP ret_data = PROTECT( allocVector( VECSXP, b_n )); 
  for(int i=0; i < b_n; ++i){
    SET_VECTOR_ELT( ret_data, i, allocVector( VECSXP, 5 ));
    SEXP rd = GET_VECTOR_ELT(ret_data, i);
    SET_VECTOR_ELT( rd, 0, allocVector( INTSXP, 1 + sizeof( struct align_stats ) / sizeof(int) ));
    SET_VECTOR_ELT( rd, 1, allocVector( STRSXP, 2 ));  // the sequences
    SET_VECTOR_ELT( rd, 2, allocVector( INTSXP, args->ret_data[i].a_pos_l ));
    SET_VECTOR_ELT( rd, 3, allocVector( INTSXP, args->ret_data[i].b_pos_l ));
    SET_VECTOR_ELT( rd, 4, allocMatrix( INTSXP, args->ret_data[i].pos_table_nrow, 2));
    // this is a bit dodgy, but trying to merge the score and the align_stats together
    INTEGER(VECTOR_ELT(rd, 0))[0] = args->ret_data[i].score;
    memcpy( (void*)(1 + INTEGER(VECTOR_ELT(rd, 0))), (void*)&(args->ret_data[i].al_stats), sizeof(struct align_stats) );
    SET_STRING_ELT( VECTOR_ELT(rd, 1), 0, mkChar( args->ret_data[i].al_a ));
    SET_STRING_ELT( VECTOR_ELT(rd, 1), 1, mkChar( args->ret_data[i].al_b ));
    memcpy( (void*)INTEGER(VECTOR_ELT(rd, 2)), args->ret_data[i].a_pos, sizeof(int) * args->ret_data[i].a_pos_l );
    memcpy( (void*)INTEGER(VECTOR_ELT(rd, 3)), args->ret_data[i].b_pos, sizeof(int) * args->ret_data[i].b_pos_l );
    memcpy( (void*)INTEGER(VECTOR_ELT(rd, 3)), args->ret_data[i].pos_table, sizeof(int) * 2 * args->ret_data[i].pos_table_nrow );
  }
  needle_wunsch_args_free( args );
  free(b_l);
  free(b_seq);
  UNPROTECT(1);
  return( ret_data );
}

static const R_CallMethodDef callMethods[] = {
  {"align_exons", (DL_FUNC)&align_exons, 8},
  {"align_seqs", (DL_FUNC)&align_seqs, 8},
  {NULL, NULL, 0}
};

void R_init_exon_aligneR(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
 
