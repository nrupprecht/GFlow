#include "Field.h"
#include "MAC.h"

int main(int argc, char* argv[]) {
  int dim = 10;
  int maxiter = 15;

  stringstream stream;
  if (argc>1) {
    stream << argv[1];
    stream >> dim;
  }

  if (argc>2) {
    stream << argv[2];
    stream >> maxiter;
  }

  MAC test(dim,dim);
  test.run(maxiter);

  return 0;

}
