//===================================================================================
// apply a prefix operation


#include "prefix.h"


//===================================================================================
// apply a prefix operator to a vector

// improvements possible in this function?

int prefix(prefix_data data, prefix_operator op, int* processor_count)
{
  *processor_count = 1;
  long index;

  for (index = 1; index < data.vector_len; index++)
    {
      char* target_pointer = data.pointer + index * op.type_len;
      op.function(target_pointer, target_pointer, target_pointer - op.type_len);
    }

  return PREFIX_SUCCESS;
}

//===================================================================================
