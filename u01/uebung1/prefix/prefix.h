//===================================================================================

// result value of prefix operation
#define PREFIX_SUCCESS (0)


// prefix operator type in C
typedef void (*prefix_operator_function)(void* r, void* a, void* b);


// operator description
typedef struct {
  char* name; // operator name
  prefix_operator_function function; // C function for operator
  void* neutral; // neutral element of operation
  int type_len; // length of one data element (in bytes)
} prefix_operator;


// data vector description
typedef struct {
  char* pointer;    // pointer to vector start
  long bytes_len;   // size of one data element (in bytes)
  long vector_len;  // length of vector (in  elements)
} prefix_data;


// interface function prototype
int prefix(prefix_data data, prefix_operator op, int* processor_count);


//===================================================================================
