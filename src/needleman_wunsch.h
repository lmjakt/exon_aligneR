#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

// in order to keep track of the starting position of
// any tracks we can put more data into the integers
// We get 15 bits of left and right which is up to
// 32000 bases; we are not going to do alignments
// that are that long. Well, we might if we had a short
// sequence against a very long one;
static const unsigned int ptr_mask =  3;
static const unsigned int left_mask = ((1 << 15) - 1) << 2;
static const unsigned int up_mask = (((1 << 15) - 1) << 2) << 15;
static const unsigned int left_shift = 2;
static const unsigned int up_shift = 17;

//int which_max_i(int *v, int l);
// try with a static inline function which_max_i
// commented to remove warning for lack of use.
/* static int which_max_i(int *v, int l){ */
/*     int max_i = 0; */
/*   double max = v[0]; */
/*   for(int i=1; i < l; i++){ */
/*     if(v[i] > max){ */
/*       max = v[i]; */
/*       max_i = i; */
/*     } */
/*   } */
/*   return(max_i); */
/* } */

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
  // sequences with gaps inserted appropriately.
  // Note that these will be 0-terminated.
  char *a_al;
  char *b_al;
  unsigned char *cigar_ops;
  int *cigar_n;
  int cigar_length;
  // pointers to alignments not masked by this alignment.
  // default to 0 values
  struct sw_alignment *top;
  struct sw_alignment *bottom;
  /* struct sw_alignment *top_left; */
  /* struct sw_alignment *top_right; */
  /* struct sw_alignment *bottom_left; */
  /* struct sw_alignment *bottom_right; */
};

// This is a 4-way recursive alignment that harvests all unique position.
// Procedure:
// 1. Find highest score within the indicated table region;
// 2. Extract the alignment and assign the values to a new sw_align
//    struct. Assign this struct to the sw_align pointer.
// 3. Call extract_sw_alignments for the four regions of the table that
//    are not shadowed by the extracted alignment. 
void extract_sw_alignments(const unsigned  char *a, const unsigned char *b,
			   int *ptr_table, int *score_table, int height, int width,
			   int row_begin, int row_end, int col_begin, int col_end,
			   struct sw_alignment **sw_align,
			   int min_width, int min_score);

void free_sw_alignments( struct sw_alignment *align );
void count_sw_alignments( struct sw_alignment *align, int *n);
void harvest_sw_aligns( struct sw_alignment *align, int *align_table,
			int **cigar_ops, int **cigar_n, int *cigar_lengths,
			int *i, int aligns_n,
			char **a_al, char **b_al);



// Aligning transcripts to genome locations with intron position
// hinting
// There are a number of problems with alignment to genomic locations:
// 1. The genomic region is too long to remember all scores and pointers.
//    This can be overcome by only remembering the previous column of scores
//    (genome sequence along the top of the matrix) and the alignments associated
//    with each score. The alignments can be encoded as cigar strings.
// 2. The score associated with horizontal moves can either indicate indels
//    or introns meaning that we need to keep track of two different scores until
//    the intron associated score has a higher value. In theory one could simply
//    have some function whereby the score goes up once the gap is long enough, but
//    one will then risk getting negative scores and giving up on the alignment.
//    Hence keeping two scores seems reasonable; any continuations to the alignment
//    before an intron length is obtained will then be considered as an insertion.
//    Continuations of the alignment after that may take into consideration the
//    expected splicing signals; as indels in theory, even within exons can be longer than
//    76 base pairs. Especially in terminal ones.

// A structure that holds a pointer
// needs to be a linked list with a reference count
// to avoid copying operations.
struct cig_op {
  // 2 smallest bits hold the cigar operation; 
  // left = 1, up = 2, diagnoal = 1 | 2 (3)
  // the three bits after that hold the child information
  // (i.e. which cells are pointing to this cell)
  unsigned char data;
  struct cig_op *parent;
};

// holds a pointer to a cig_op chain and some information about that chain
// to allow the tracing of an alignment
struct cig_data {
  struct cig_op *ops;
  int score;
  // current position (i.e. the end of the operations)
  int pos_t;
  int pos_g;
  // positions related to the alignment
  int beg_t;
  int beg_g;
  // positions related to the peak;
  int max_score;
  int max_t;
  int max_g;
  // the number of consecutive left ops
  int left_n;
  int up_n;
};

// the constant mallocing of data is likely to be slow
// as we have to do this every time we extend or start
// a new alignment (and this willl be very frequently)
// I have some ideas as to how I can do something more
// efficiently, but we should try the straightforward way
// first..
struct cig_op *init_cig_op(unsigned char op, struct cig_op *parent);

unsigned char c_op(struct cig_op *co);

void clear_ops(struct cig_op *node, unsigned char op);

// make a vector of cig_data structures
// with all values set to 0
struct cig_data *init_cig_data();

// A matrix of DP pointers encoded as 2 bits in a byte
// row minor (as in R)
// ptr[0] is stored in the least signficant bit of a byte
// note that a ptr can only be set once! This is to make the
// the set_ptr operation faster. This may be changed in the future.
struct bit_ptr {
  unsigned char *data;
  unsigned int nrow;
  unsigned int ncol;
};

struct bit_ptr init_bit_ptr(unsigned int nrow, unsigned int ncol);
void free_bit_ptr(struct bit_ptr *ptr);
unsigned char get_ptr(struct bit_ptr *ptr, unsigned int row, unsigned int col);
void set_ptr(struct bit_ptr *ptr, unsigned int row, unsigned int col, unsigned char val);

// sw_alignment is actually intended to hold a linked list of
// sw_alignment; however here we will use it to hold a single alignment
// extractd from cig_data;
// returns 0 on success
int extract_swi_alignment( struct cig_data *cd, const unsigned char *tr, const unsigned char *gen,
			   struct sw_alignment *alignment);

int extract_bitp_alignment( const unsigned char *tr, const unsigned char *gen, 
			    struct bit_ptr *ptr, int max_row, int max_col, int score,
			    struct sw_alignment *alignment);

// I haven't yet decided what this function should return.
// It is expected that tr is a sequence of characters, one of which
// represents introns (by default I). This should not be alignable to any
// character but must be passed by gaps in the tr sequence (if intron is
// specific to the transcript) or in the gen sequence (if they are shared).
// In addition, introns specific to the genome may also exist; possibly
// we should have some form of penalty for such introns, similar to ones
// specific to the transcript.

// The innovation in this function is to evaluate the impact of gaps on the score
// at the point when the gap is joined; i.e. the gap is free until we try to align
// characters again at which point we consider how to penalise the gap
struct sw_alignment *align_transcript_to_genome( const unsigned char *tr, const unsigned char *gen,
						int tr_l, int gen_l, int gap_i, int gap_e, int intron_p,
						int *sub_table, int al_offset, int al_size, char intron_char );


// a function to align making use of the bit_ptr
// initially make it void; later on work out what return values to provide

struct sw_alignment *align_transcript_bitp(const unsigned char *tr, const unsigned char *gen,
					   int tr_l, int gen_l, int gap_i, int gap_e, int intron_p,
					   int *sub_table, int al_offset, int al_size, char intron_char);


struct sw_alignment *align_transcript_bitp_ssa(const unsigned char *tr, const unsigned char *gen,
					       int tr_l, int gen_l, float gap_i, float gap_e_mod, float fshift,
					       float intron_new, float intron_bad, int min_intron_size,
					       int *sub_table, int al_offset, int al_size, char intron_char,
					       int global);

static inline float gap_penalty(unsigned int gap_l, float min_is, float gap_i, float mod, float fshift[3]);

// A struct to hold the set of penalties that we wish to use:
struct al_penalty {
  int *sub_table;
  int al_offset;
  int al_size;
  char intron_char;
  float gap_i;
  // gap_e_mod modifies the extension: the bigger the cheaper the extension as it divides
  float gap_e_mod; 
  float fshift[3];
  float intron_new;
  float intron_bad;
  int min_intron_size;
  // intron starts and ends:
  int intron_type_n;
  // precalculated cumulative gap values for 0 -> min_intron_size;
  // (size of array = min_intron_size + 1);
  float *gap;
  const char **intron_start;
  const char **intron_end;
};

// this is just a holder; it should not allocate or free anything
struct al_penalty init_al_penalty(int *sub_table, int al_offset, int al_size, char intron_char,
				  float gap_i, float gap_e_mod, float fshift, float intron_new, float intron_bad, int min_intron_size,
				  int intron_type_n, const char **intron_start, const char **intron_end
				  );

void free_al_penalty(struct al_penalty *pen);

#endif
