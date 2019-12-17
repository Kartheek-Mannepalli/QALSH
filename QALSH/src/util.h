#ifndef __UTIL_H
#define __UTIL_H


// -----------------------------------------------------------------------------
//  Global variables
// -----------------------------------------------------------------------------
extern long g_memory;

// -----------------------------------------------------------------------------
//  Uitlity functions
// -----------------------------------------------------------------------------
int compfloats(						// compare two float values
	float v1,							// 1st float value
	float v2);							// 2nd float value

// -----------------------------------------------------------------------------
//void error(							// an error message
//	char* msg,							// an message
//	bool is_exit);						// whether exit the program

// -----------------------------------------------------------------------------
int check_mem();					// check memory is enough

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L2 distance (type float)
	float* p1,							// 1st point
	float* p2,							// 2nd point
	int dim);							// dimension


// -----------------------------------------------------------------------------
//  Functions used for the input/output of datasets and query sets.
// -----------------------------------------------------------------------------
int read_set(						// read (data or query) set from disk
	int n,								// number of data points
	int d,								// dimensionality
	char* set,							// address of dataset
	float** points);					// data or queries (return)

// -----------------------------------------------------------------------------
int write_data_new_form(			// write dataset with new format
	int n,								// cardinality
	int d,								// dimensionality
	int B,								// page size
	float** data,						// data set
	char* output_path);					// output path

// -----------------------------------------------------------------------------
void get_data_filename(				// get file name of data
	int data_id,						// data file id
	char* data_path,					// path to store data in new format
	char* fname);						// file name of data (return)

// -----------------------------------------------------------------------------
void write_data_to_buffer(			// write data to buffer
	int d,								// dimensionality
	int left,							// left data id
	int right,							// right data id
	float** data,						// data set
	char* buffer);						// buffer to store data

// -----------------------------------------------------------------------------
int write_buffer_to_page(			// write data to one page
	int B,								// page size
	char* fname,						// file name of data
	char* buffer);						// buffer to store data

// -----------------------------------------------------------------------------
int read_data(						// read data from page
	int id,								// index of data
	int d,								// dimensionality
	int B,								// page size
	float* data,						// real data (return)
	char* output_path);					// output path

// -----------------------------------------------------------------------------
int read_buffer_from_page(			// read data from page
	int B,								// page size
	char* fname,						// file name of data
	char* buffer);						// buffer to store data

// -----------------------------------------------------------------------------
void read_data_from_buffer(			// read data from buffer
	int index,							// index of data in buffer
	int d,								// dimensionality
	float* data,						// data set
	char* buffer);						// buffer to store data

void get_leading_folder(char *_path, char *_folder); // omid

#endif

// Detailed time and IO analysis (omid) (begin)
#ifndef GLOBALVARS_INCLUDED
#define GLOBALVARS_INCLUDED
#ifndef GLOBALVARS_CPP
#define DECLARE_EXTERN extern
#define EXTERN_INIT(x)
#else
#define DECLARE_EXTERN
#define EXTERN_INIT(x) x
#endif
DECLARE_EXTERN double INDEXING_TIME;
DECLARE_EXTERN double GT_TIME;
DECLARE_EXTERN double QUERY_TIME;
DECLARE_EXTERN double ALG_TIME;
DECLARE_EXTERN double DIST_CALC_TIME;
DECLARE_EXTERN double INDEX_IO_TIME;
DECLARE_EXTERN double DATA_IO_TIME;
DECLARE_EXTERN int INDEX_IO_NUM EXTERN_INIT(= 0);
DECLARE_EXTERN int DATA_IO_NUM EXTERN_INIT(= 0);
DECLARE_EXTERN double INDEX_IO_SIZE;
DECLARE_EXTERN double DATA_IO_SIZE;
DECLARE_EXTERN char INDEX_PATH[1000];
DECLARE_EXTERN char DATA_BIN_PATH[1000];
DECLARE_EXTERN bool isQinDS;
DECLARE_EXTERN double VR_ALG_TIME;
DECLARE_EXTERN double VR_IO_TIME;
DECLARE_EXTERN double VR_TIME;
DECLARE_EXTERN int VR_COUNT EXTERN_INIT(= 0);
DECLARE_EXTERN bool VR_FLAG EXTERN_INIT(= false);

// Timers for indexing (omid) (begin)
DECLARE_EXTERN double READ_DS_TIME;
DECLARE_EXTERN double WRITE_DS_BIN_TIME;
DECLARE_EXTERN double INIT_PARAMS_TIME;
DECLARE_EXTERN double INIT_HASH_TIME;
DECLARE_EXTERN double PROJ_POINTS_TIME;
DECLARE_EXTERN double SORT_TIME;
DECLARE_EXTERN double BUILD_WRITE_TREE_TIME;
// Timers for indexing (omid) (end)


#endif

// Detailed time and IO analysis (omid) (end)
