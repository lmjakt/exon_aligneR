#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

int which_max_i(int *v, int l);
int which_max_d(double *v, int l);
int m_offset(int row, int column, int height);

void needleman_wunsch( const unsigned char *a, const unsigned char *b, int a_l, int b_l,
		       int gap_i, int gap_e,
		       int *sub_table, int al_offset, int al_size,
		       int tgaps_free,
		       int *score_table, int *ptr_table);

void extract_nm_alignment(int* pointers, int height, int width, const char *a, const char *b,
			  char **a_a, char **b_a);

void char_at(const char *word, char c, int **pos, int *pos_l);
int *aligned_i(int *pos1, int *pos2, int l1, int l2, int *nrow);

#endif
