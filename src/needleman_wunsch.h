#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

//int which_max_i(int *v, int l);
// try with a static inline function which_max_i
static int which_max_i(int *v, int l){
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

int which_max_d(double *v, int l);
int m_offset(int row, int column, int height);

struct align_stats {
  int al_length;  // the total alignment including gaps
  int al_n;       // the number of aligned positions in a
  int a_gap_i;    // number of gap insertions in the first sequence
  int b_gap_i;    // number of gap insertions in the second sequence
  int a_gap;      // total number of gaps in sequence
  int b_gap;      //
  int a_gap_l;    // terminal gaps
  int a_gap_r;    // terminal gaps
  int b_gap_l;    // terminal gaps
  int b_gap_r;    // terminal gaps
  int match_n;
  int mismatch_n;
  int transition;
  int transversion;
  int A, C, G, T;  // individual counts..
};

struct align_stats align_stats_init();

void needleman_wunsch( const unsigned char *a, const unsigned char *b, int a_l, int b_l,
		       int gap_i, int gap_e,
		       int *sub_table, int al_offset, int al_size,
		       int tgaps_free,
		       int *score_table, int *ptr_table);

struct align_stats extract_nm_alignment(int* pointers, int height, int width, const unsigned char *a, const unsigned char *b,
			  char **a_a, char **b_a);

void char_at(const char *word, char c, int **pos, int *pos_l);
int *aligned_i(int *pos1, int *pos2, int l1, int l2, int *nrow);

// do we have a transition or a transversion.
// -1 => transition
//  0 => unknown
// +1 = transversion
// a and b must be in A,C,G,T to give non-0 result
int mut_type(char a, char b, int *A, int *C, int *G, int *T);

// Functions for local alignment. We should either rename this file to
// something like, dynamic_alignment.h, or create a separate
// header file for smith_waterman, and a header for common data structures
// But one step at a time..

// smith_waterman should probably return something that tells us the location
// of the maximall scoring cell; but let us not care too much

void smith_waterman( const unsigned char *a, const unsigned char *b, int a_l, int b_l,
		     int gap_i, int gap_e,
		     int *sub_table, int al_offset, int al_size,
		     int *score_table, int *ptr_table,
		     int *max_score, int *max_row, int *max_column);


// A linked list that can be used when extracting all non-masked alignments
// From a given region.

struct sw_alignment {
  int row_begin, row_end;
  int col_begin, col_end;
  int score;
  // the length of the alignment including inserts and deleteions
  int al_length;
  // We encode the cigar string with one char and one integer per
  // operation (in two arrays). This is incredibly wasteful, but this
  // can be directly handled in R.
  unsigned char *cigar_ops;
  int *cigar_n;
  int cigar_length;
  // pointers to alignments not masked by this alignment.
  // default to 0 values
  struct sw_alignment *top_left;
  struct sw_alignment *top_right;
  struct sw_alignment *bottom_left;
  struct sw_alignment *bottom_right;
};

// This is a 4-way recursive alignment that harvests all unique position.
// Procedure:
// 1. Find highest score within the indicated table region;
// 2. Extract the alignment and assign the values to a new sw_align
//    struct. Assign this struct to the sw_align pointer.
// 3. Call extract_sw_alignments for the four regions of the table that
//    are not shadowed by the extracted alignment. 
void extract_sw_alignments(int *ptr_table, int *score_table, int height, int width,
			   int row_begin, int row_end, int col_begin, int col_end,
			   struct sw_alignment **sw_align);

void free_sw_alignments( struct sw_alignment *align );
void count_sw_alignments( struct sw_alignment *align, int *n);

#endif
