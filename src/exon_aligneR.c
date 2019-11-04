#include <R.h>
#include <Rinternals.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>


// This is a rewrite of exon_aligneR.c
// Modified to:
// 1. Use filly by column matrices like in R
// 2. Always place seq 1 (a) by rows and seq 2 (b) by columns
// 3. Remove the merge functionaliy

// To enforce the matrices; always use this
int m_offset(int row, int column, int height){
  return( column * height + row );
}


// Note that we no longer use the exon_set here. We instead use a transcript
// As that makes it much easier to merge adjoining exons.
struct exon_set {
  int n;
  char **seq;
  int *seq_l;
};

// seq_r must be STRSXP of non-zero length
// This makes a copy of the strings. That ought not to be necessary, but... 
struct exon_set exon_set_init(SEXP seq){
  struct exon_set e_set;
  e_set.n = length(seq);
  e_set.seq = malloc( sizeof(const char*) * e_set.n );
  e_set.seq_l = malloc( sizeof(int) * e_set.n );
  for(int i=0; i < e_set.n; ++i){
    SEXP s = STRING_ELT(seq, i);
    e_set.seq_l[i] = length(s);
    e_set.seq[i] = malloc( sizeof(char) * (length(s) + 1) );
    strncpy( e_set.seq[i], CHAR(s), length(s) + 1);
    e_set.seq[i][ length(s) ] = 0;
    //    e_set.seq[i] = CHAR(s);
  }
  return(e_set);
}

void exon_set_free(struct exon_set e_set){
  for(int i=0; i < e_set.n; ++i)
    free( e_set.seq[i] );
  free( e_set.seq );
  free( e_set.seq_l );
}


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


int which_max(double *v, int l){
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


// takes pre-allocated tables
void init_nm_tables( double *scores, char *pointers, int m_height, int m_width,
		     double gap_i, double gap_e, const char left, const char up){
  memset( (void*)scores, 0, sizeof(double) * m_height * m_width );
  memset( (void*)pointers, 0, sizeof(char) * m_height * m_width);

  // Set the initial scores and pointers
  scores[ m_offset(0, 1, m_height) ] = gap_i;
  scores[ m_offset(1, 0, m_height) ] = gap_i;
  pointers[ m_offset(0, 1, m_height) ] = left;
  pointers[ m_offset(1, 0, m_height) ] = up;

  // and then fill in the first row and column
  for(int i=2; i < m_width; ++i){
    scores[ m_offset(0, i, m_height) ] = scores[ m_offset(0, i-1, m_height) ] + gap_e;
    pointers[ m_offset(0, i, m_height) ] = left;
  }
  for(int i=2; i < m_height; ++i){
    scores[ m_offset(i, 0, m_height) ] = scores[ m_offset(i-1, 0, m_height) ] + gap_e;
    pointers[ m_offset(i, 0, m_height) ] = up;
  }

}

// expects to be given a table for the pointers since this is necessary to trace the alignment back
double exon_nm(const char *seq_a, const char *seq_b, int a_l, int b_l, double *penalties, double *scores, char *pointers){
  //double exon_nm(struct exon_set a, struct exon_set b, int a_i, int b_i, double *penalties, double *scores, char *pointers){
  /* if(a_i >= a.n || b_i >= b.n) */
  /*   error("a_i >= a.n or b_i >= b.n"); */
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
  init_nm_tables( scores, pointers, m_height, m_width, gap_i, gap_e, left, up);
  
  // And then we simply go through the table positions. Let us make a pointer to the two
  // sequences that we are using
  /* char *seq_a = a.seq[a_i]; */
  /* char *seq_b = b.seq[b_i]; */

  int o=0, o_l, o_u, o_d;  // offsets for the different positions
  
  for(int row=1; row < m_height; ++row){
    for(int column=1; column < m_width; ++column){
      double sc[3];
      o = m_offset(row, column, m_height);
      o_l = m_offset(row, column-1, m_height);
      o_u = m_offset(row-1, column, m_height);
      o_d = m_offset(row-1, column-1, m_height);
      sc[0] = scores[o_l] + ( pointers[o_l] == left ? gap_e : gap_i );
      sc[1] = scores[o_u] + ( pointers[o_u] == up ? gap_e : gap_i );
      sc[2] = scores[o_d] + (seq_a[row-1] == seq_b[column-1] ? match : mis_match );
      int max_i = which_max(sc, 3);
      scores[o] = sc[max_i];
      pointers[o] = max_i + 1;
    }
  }
  return(scores[o]);
}

// the use of int for the pointers is wasteful, but it makes it easier to
// interface with R. I would otherwise have to build a table of char*,
// which is even worse.. 
/* void gene_nm( struct exon_set a, struct exon_set b, double *exon_scores, */
/* 	      double match, double gap, */
/* 	      double *scores, int *pointers ){ */
void gene_nm( struct transcript a, struct transcript b, double *exon_scores,
	      double match, double gap,
	      double *scores, int *pointers ){
  int height = a.e_n + 1;
  int width = b.e_n + 1;
  int m_size = width * height;
  memset( (void*)scores, 0, sizeof(double) * m_size );
  memset( (void*)pointers, 0, sizeof(int) * m_size );

  int left = 1;
  int up = 2;
  // init the pointers..
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

  for(int row=1; row < height; ++row){
    for(int column=1; column < width; ++column){
      int o = m_offset(row, column, height);
      int o_l = m_offset(row, column-1, height);
      int o_u = m_offset(row-1, column, height);
      int o_d = m_offset(row-1, column-1, height);
      double sc[3];
      sc[0] = scores[o_l] - match * b.e_lengths[column-1] * gap;
      sc[1] = scores[o_u] - match * a.e_lengths[row-1] * gap;
      sc[2] = scores[o_d] + exon_scores[ m_offset(row-1, column-1, height-1) ];
      int max_i = which_max( sc, 3 );
      scores[o] = sc[ max_i ];
      pointers[o] = max_i + 1;
    }
  }
}
	      

// a_seq_r: the sequences of exons of gene a
// b_seq_r: the sequences of exons of gene b
// e_penalties_r: the penalties / scores used for aligning sequences to each other
// g_penalties_r: some parameters that can be used to define the penalties
//                for the final exon alignments
//SEXP align_exons(SEXP a_seq_r, SEXP b_seq_r, SEXP e_penalties_r, SEXP g_penalty_r ){
SEXP align_exons(SEXP a_seq_r, SEXP b_seq_r, SEXP a_lengths, SEXP b_lengths,
		 SEXP e_penalties_r, SEXP g_penalty_r ){
  // SANITY check
  /* if( TYPEOF(a_seq_r) != STRSXP || TYPEOF(b_seq_r) != STRSXP ) */
  /*   error("The first two arguments should be character vectors"); */
  if( TYPEOF(e_penalties_r) != REALSXP || TYPEOF(g_penalty_r) != REALSXP )
    error("Arguments 3 and 4 should vectors of real numbers");
  if(length(e_penalties_r) != 4)
    error("The third argument should provide: match, mismatch, gap insertion, gap extension values");
  if(length(g_penalty_r) != 1)
    error("The fourth argument should provide a single value from which we can derive an exon gap penalty");
  
  struct transcript trans_a = transcript_init(a_seq_r, a_lengths);
  struct transcript trans_b = transcript_init(b_seq_r, b_lengths);
  double *e_penalties = REAL(e_penalties_r);
  double g_penalty = REAL(g_penalty_r)[0];

  /* // let us also define the number of exons */
  /* int a_n = length(a_seq_r); */
  /* int b_n = length(b_seq_r); */
  /* if(a_n < 1 || b_n < 1) */
  /*   error("There must be at least one exon for each sequence"); */

  /* // Define two exon sets to hold the data */
  /* struct exon_set a_exons = exon_set_init( a_seq_r ); */
  /* struct exon_set b_exons = exon_set_init( b_seq_r ); */

  // Allocate space for a table of scores that we can return to R
  int a_n = trans_a.e_n;
  int b_n = trans_b.e_n;
  SEXP ret_data = PROTECT( allocVector( VECSXP, 3 ));
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
      exon_scores[ m_offset(row, column, a_n) ] = exon_nm( trans_a.seq + trans_a.e_offsets[row],
							   trans_b.seq + trans_b.e_offsets[column],
							   trans_a.e_lengths[row], trans_b.e_lengths[column],
							   e_penalties, scores, pointers );
      /* exon_scores[ m_offset(row, column, a_n) ] = exon_nm( a_exons.seq[row], b_exons.seq[column], */
      /* 							   a_exons.seq_l[row], b_exons.seq_l[column], */
      /* 							   e_penalties, scores, pointers ); */
    }
  }

  //  gene_nm( a_exons, b_exons, exon_scores, e_penalties[0], g_penalty, gene_score_matrix, gene_pointer_matrix );
  gene_nm( trans_a, trans_b, exon_scores, e_penalties[0], g_penalty, gene_score_matrix, gene_pointer_matrix );
  
  free(scores);
  free(pointers);
  transcript_free( trans_a );
  transcript_free( trans_b );
  /* exon_set_free( a_exons ); */
  /* exon_set_free( b_exons ); */
  UNPROTECT(1);
  return(ret_data);
}
