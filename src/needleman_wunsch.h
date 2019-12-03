#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

int which_max_i(int *v, int l);
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
};

struct align_stats align_stats_init();

void needleman_wunsch( const unsigned char *a, const unsigned char *b, int a_l, int b_l,
		       int gap_i, int gap_e,
		       int *sub_table, int al_offset, int al_size,
		       int tgaps_free,
		       int *score_table, int *ptr_table);

struct align_stats extract_nm_alignment(int* pointers, int height, int width, const char *a, const char *b,
			  char **a_a, char **b_a);

void char_at(const char *word, char c, int **pos, int *pos_l);
int *aligned_i(int *pos1, int *pos2, int l1, int l2, int *nrow);


#endif
