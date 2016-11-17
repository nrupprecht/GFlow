#include "GField.h"

Index::Index(int first) {
  total = 1;
  entries = new int(first);
}

GField::GField() : spacing(0), stride(0), array(0), lowerBounds(0), upperBounds(0) {};

GField::~GField() {
  /*
  if (spacing) delete [] spacing;
  if (stride) delete [] stride;
  if (array) delete [] array;
  if (lowerBounds) delete [] lowerBounds;
  if (upperBounds) delete [] upperBounds;
  */
}

void GField::initialize(Shape s) {
  shape = s;
  total = s.getTotal();
  if (total==0) return;
  // Set stride array
  int count = 1, rank = shape.rank;
  stride = new int[rank];
  for (int i=0; i<rank; i++) {
    count *= shape.dims[i];
    stride[i] = total/count;
  }
}

GField& GField::operator+=(const GField& F) {
  if (F.shape!=shape) throw GFieldMismatch();
  for (int i=0; i<total; i++)
    array[i] += F.array[i];
}
