#include <R.h>
#include <Rinternals.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>

// Multithreading to be implemented at a later stage;
// First work out what functions will be needed to implement
// The algorithm.

struct exon_set {
  int n;
  const char **seq;
  int *seq_l;
};

// seq_r must be STRSXP of non-zero length
struct exon_set exon_set_init(SEXP seq){
  struct exon_set e_set;
  e_set.n = length(seq);
  e_set.seq = malloc( sizeof(const char*) * e_set.n );
  e_set.seq_l = malloc( sizeof(int) * e_set.n );
  for(int i=0; i < e_set.n; ++i){
    SEXP s = STRING_ELT(seq, i);
    e_set.seq_l[i] = length(s);
    e_set.seq[i] = CHAR(s);
  }
  return(e_set);
}

void exon_set_free(struct exon_set e_set){
  free( e_set.seq );
  free( e_set.seq_l );
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

// align the sequence starting at exon a_i for a_n exons from a
// with the sequence staring at exon b_i for b_n exons from b
// This function is written in order to allow the fusion of exons
// So that I can align a single exon to two separate exons in a different sequence
// For this purpose either a_n or b_n should always be 1;
// penalties must be:
// match, mismatch, gap insertion, gap extension
double exon_nm(struct exon_set a, struct exon_set b, int a_i, int b_i, int a_n, int b_n,
	       double *penalties){
  // for ease of reading and writing we do:
  double match = penalties[0];
  double mis_match = penalties[1];
  double gap_i = penalties[2];
  double gap_e = penalties[3];
  // the following are for the pointers
  const char left = 1;
  const char up = 2;
  //  const char diag = 3;
  int a_l = 0;
  int b_l = 0;
  for(int i=a_i; i < a_i + a_n; ++i)
    a_l += a.seq_l[i];
  for(int i=b_i; i < b_i + b_n; ++i)
    b_l += b.seq_l[i];
  // I then have to keep careful track of the position in the sequences.
  // Since we only want the final score we only need two rows of data where
  // we update one after the other...
  double *all_scores = malloc(sizeof(double) * (a_l + 1) * 2);
  memset( (void*)all_scores, 0, sizeof(double) * (a_l + 1) * 2);
  double *scores[2] = { all_scores, all_scores + a_l + 1 };
  // Since we want to use affine gap penalties, we will also need
  // to keep track of the pointers. We can use a char array in the same way.
  char *all_ptr = malloc(sizeof(char) * (a_l + 1) * 2);
  memset( (void*)all_ptr, 0, sizeof(char) * (a_l + 1) * 2);
  char *ptr[2] = {all_ptr, all_ptr + a_l + 1};
  // set the scores of the first set of rows. 
  scores[0][1] = gap_i;
  for(int i=2; i <= a_l; ++i){
    scores[0][i] = scores[0][i-1] + gap_e;
    ptr[0][i] = left;
  }
  // iterate through the row s
  // we need the local offsets for the sequence positions
  int a_o = 0;
  int a_j = 0;
  int b_o = 0;
  int r_1 = 0;
  int r_2 = 1;
  for(int row=1; row <= b_l; row++){
    if(b_o == b.seq_l[b_i]){
      b_o = 0;
      b_i++;
    }
    char b_res = b.seq[b_i][b_o];
    r_1 = (1 + row) % 2;
    r_2 = row % 2;
    // r_2 is the row that we are working on. We first need to set the first value of r_2
    scores[r_2][0] = (row == 1) ? gap_i : scores[r_1][0] + gap_e;
    ptr[r_2][0] = up;
    for(int column=1; column <= a_l; column++){
      if(a_o == a.seq_l[a_j]){
	a_o = 0;
	a_j++;
      }
      char a_res = a.seq[a_j][a_o];
      // and here we do the actual align.
      double cell_scores[3];
      cell_scores[0] = (ptr[r_2][column - 1] == left) ? scores[r_2][column - 1] + gap_e : scores[r_2][column - 1] + gap_i;
      cell_scores[1] = (ptr[r_1][column] == up) ? scores[r_1][column] + gap_e : scores[r_1][column] + gap_i;
      cell_scores[2] = (a_res == b_res) ? scores[r_1][column-1] + match : scores[r_1][column-1] + mis_match;
      int max_i = which_max(cell_scores, 3);
      scores[r_2][column] = cell_scores[max_i];
      ptr[r_2][column] = max_i + 1;
      a_o++;
    }
    a_o = 0;
    a_j = a_i;
    b_o++;
  }
  double max_score = scores[r_2][a_l];
  free(all_scores);
  free(all_ptr);
  return(max_score);
}

// a_seq_r: the sequences of exons of gene a
// b_seq_r: the sequences of exons of gene b
// e_penalties_r: the penalties / scores used for aligning sequences to each other
// g_penalties_r: some parameters that can be used to define the penalties
//                for the final exon alignments
SEXP align_exons(SEXP a_seq_r, SEXP b_seq_r, SEXP e_penalties_r, SEXP g_penalty_r ){
  // First, define the alignment scores for all exons vs all exons in order to
  // then be able to use:
  // log( score / mean(scores) )
  // or possibly
  // log( score_{i,j} / mean( [score_{i}, score_{j}] ) )
  // as the score / penalty for aligning a given exon-exon pair
  // Then create an exon alignment considering each exon as a single residue

  // SANITY check
  if( TYPEOF(a_seq_r) != STRSXP || TYPEOF(b_seq_r) != STRSXP )
    error("The first two arguments should be character vectors");
  if( TYPEOF(e_penalties_r) != REALSXP || TYPEOF(g_penalty_r) != REALSXP )
    error("Arguments 3 and 4 should vectors of real numbers");
  if(length(e_penalties_r) != 4)
    error("The third argument should provide: match, mismatch, gap insertion, gap extension values");
  if(length(g_penalty_r) != 1)
    error("The fourth argument should provide a single value from which we can derive an exon gap penalty");
  
  double *e_penalties = REAL(e_penalties_r);
  double g_penalty = REAL(g_penalty_r)[0];

  // let us also define the number of exons
  int a_n = length(a_seq_r);
  int b_n = length(b_seq_r);
  if(a_n < 1 || b_n < 1)
    error("There must be at least one exon for each sequence");

  // Define two exon sets to hold the data
  struct exon_set a_exons = exon_set_init( a_seq_r );
  struct exon_set b_exons = exon_set_init( b_seq_r );

  // We then wish to make a matrix alignment scores for each exon vs each exon...
  // This could easily be done in parallell
  SEXP exon_al_scores_r = PROTECT(allocMatrix( REALSXP, a_n, b_n ));
  double *exon_al_scores = REAL( exon_al_scores_r );
  //  double *exon_al_scores = malloc( sizeof(double) * a_n * b_n );
  for(int i=0; i < a_n; ++i){
    for(int j=0; j < b_n; ++j){
      exon_al_scores[ i * b_n + j ] = exon_nm( a_exons, b_exons, i, j, 1, 1, e_penalties );
    }
  }
  // At this point to test the the methods we could simply return a matrix containing the e_penalties
  UNPROTECT(1);
  return( exon_al_scores_r );
}
