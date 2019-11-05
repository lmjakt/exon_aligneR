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

  /* // We can use this for local alignments by setting gap_i and gap_e to 0. If that is the */
  /* // the case it is more reasonable to return from here. */
  /* if(gap_i == 0 && gap_e == 0) */
  /*   return; */
  
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
  if(!local)
    init_nm_tables( scores, pointers, m_height, m_width, gap_i, gap_e, left, up);
  else
    init_nm_tables( scores, pointers, m_height, m_width, 0, 0, left, up);
  
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
      // To 
      left_gap = (local && row == m_height - 1) ? 0 : ( pointers[o_l] == left ? gap_e : gap_i );
      up_gap = (local && column == m_width - 1) ? 0 : ( pointers[o_u] == up ? gap_e : gap_i );
      
      sc[0] = scores[o_l] + left_gap; // ( pointers[o_l] == left ? gap_e : gap_i );
      sc[1] = scores[o_u] + up_gap; // ( pointers[o_u] == up ? gap_e : gap_i );
      sc[2] = scores[o_d] + (seq_a[row-1] == seq_b[column-1] ? match : mis_match );
      int max_i = which_max(sc, 3);
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
    //  }
  
  // if the alignment is local we also need to determine the maximum bottom
  /* for(int column=1; column < m_width; ++column){ */
  /*   o = m_offset(m_height-1, column, m_height); */
  /*   dp_max_update( &max, scores[o], m_height - 1, column ); */
  /* } */
  /* return( max ); */
}

// the use of int for the pointers is wasteful, but it makes it easier to
// interface with R. I would otherwise have to build a table of char*,
// which is even worse.. 
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

// height and width are the dimensions of the table, not the number of exons
void extract_gene_alignment(int *pointers, int height, int width, struct gene_alignment* align){
  int al_length = 0;
  int row = height -1;
  int column = width -1;
  while(row > 0 || column > 0){
    int o = m_offset( row, column, height );
    if(pointers[o] == 0)
      break;
    row = (pointers[o] &  2) ? row - 1 : row;
    column = (pointers[o] & 1) ? column - 1 : column;
    al_length++;
  }
  align->length = al_length;
  align->alignment = malloc( sizeof(int) * 2 * al_length );
  int *al_a = align->alignment;
  int *al_b = align->alignment + al_length;
  row = height -1;
  column = width -1;
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
		 SEXP e_penalties_r, SEXP g_penalty_r, SEXP local_r ){
  // SANITY check
  if( TYPEOF(e_penalties_r) != REALSXP || TYPEOF(g_penalty_r) != REALSXP )
    error("Arguments 3 and 4 should vectors of real numbers");
  if(length(e_penalties_r) != 4)
    error("The third argument should provide: match, mismatch, gap insertion, gap extension values");
  if(length(g_penalty_r) != 1)
    error("The fourth argument should provide a single value from which we can derive an exon gap penalty");
  if(TYPEOF(local_r) != LGLSXP)
    error("The fifth argument should be a logical value determining whether we use localised NM-alignments");
  
  // the transcript_init function sanity checks its SEXP arguments
  struct transcript trans_a = transcript_init(a_seq_r, a_lengths);
  struct transcript trans_b = transcript_init(b_seq_r, b_lengths);
  double *e_penalties = REAL(e_penalties_r);
  double g_penalty = REAL(g_penalty_r)[0];
  char local = (char)asLogical(local_r);
  
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

  gene_nm( trans_a, trans_b, exon_scores, e_penalties[0], g_penalty, gene_score_matrix, gene_pointer_matrix );
  struct gene_alignment g_align;
  extract_gene_alignment( gene_pointer_matrix, trans_a.e_n + 1, trans_b.e_n + 1, &g_align );
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
