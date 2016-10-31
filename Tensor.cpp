#include "Tensor.h"

Tensor::Tensor(Shape s) {
  initialize(s);
}

Tensor::Tensor(const Tensor& T) {
  initialize(T.shape);
  for (int i=0; i<total; i++) array[i] = T.array[i];
}

// Assumes all the tensors have the same shape
Tensor::Tensor(const vector<Tensor>& Tvec) {
  int first = Tvec.size();
  if (first==0) {
    initialize(Shape());
    return;
  }
  Shape s = Shape(first, Tvec.at(0).getShape());
  initialize(s);
  int i=0, tot=Tvec.at(0).getShape().getTotal();
  for (auto T : Tvec) {
    for (int j=0; j<tot; j++) array[i*tot+j] = T.array[j];
    i++;
  }
}

Tensor::Tensor(const vector<double>& vec) {
  Shape s(vec.size());
  initialize(s);
  for (int i=0; i<vec.size(); i++) array[i] = vec.at(i);
}

Tensor::~Tensor() {
  if (array) delete [] array;
  if (stride) delete [] stride;
}

Tensor& Tensor::operator=(const Tensor& T) {
  if (T.total!=total && array) delete [] array; // Need a new array
  if (stride) delete [] stride;
  
  initialize(T.shape);
  // Set values
  for (int i=0; i<total; i++) array[i] = T.array[i];
}

/*
Tensor& Tensor::operator=(const Tensor&& T) {
  if (array) delete [] array;
  if (stride) delete [] stride;
  array = T.array;
  stride = T.stride;
  shape = T.shape;
  total = T.total;
}
*/

Tensor Tensor::shift(const Shape& shift) const {
  Tensor T(shape); // Same shape, shift entries
  vector<int> point;
  point.reserve(shape.rank);
  shift_helper(shift, point, T);
  return T;
}

void Tensor::set(double value, vector<int> indices, const Shape& shift) {
  if (indices.size()!=shape.rank || shift.rank!=shape.rank) return;
  for (int i=0; i<indices.size(); i++) {
    int& add = indices.at(i) = indices.at(i)+shift.at(i);
    if (add<0 || add>=shape.dims[i]) return;
  }
  at(indices) = value;
}

void multiply(const Tensor& A, int aI, const Tensor& B, int bI, Tensor& C) {
  // Check index correctness
  if (aI>=A.shape.rank || bI>=B.shape.rank) throw Tensor::TensorBadContraction();
  // Check Matrix compatability
  if (A.shape.rank + B.shape.rank - 2 != C.shape.rank) throw Tensor::TensorRankMismatch();

  int i, j;
  for (i=0, j=0; i<A.shape.rank; i++) {
    if (i==aI) continue;
    if (A.shape.dims[i]!=C.shape.dims[j]) throw Tensor::TensorDimsMismatch();
    j++;
  }
  for (i=0; i<B.shape.rank; i++) {
    if (i==bI) continue;
    if (B.shape.dims[i]!=C.shape.dims[j]) throw Tensor::TensorDimsMismatch();
    j++;
  }

  // Special case - (m,k) times (k,n) -> (m, n)
  if (A.shape.rank==2 && B.shape.rank==2) { // Normal matrix multiplication
    int m, n, k, ac, bc, cc;

    auto AT = CblasTrans;
    if (aI==0) { // A is "transposed"
      m = A.shape.dims[1];
      k = A.shape.dims[0];
    }
    else {
      m = A.shape.dims[0];
      k = A.shape.dims[1];
      AT = CblasNoTrans;
    }
    ac = A.shape.dims[1];

    auto BT = CblasNoTrans;
    if (bI==0) n = B.shape.dims[1]; // B not "transposed"
    else {
      n = B.shape.dims[0];
      BT = CblasTrans;
    }
    bc = B.shape.dims[1];
    cc = C.shape.dims[1];
    double ALPHA = 1.0, BETA = 0;

    cblas_dgemm(CblasRowMajor, AT, BT, m, n, k, ALPHA, A.array, ac, B.array, bc, BETA, C.array, cc);
    return;
  }
  if (A.shape.rank==2 && B.shape.rank==1) { // Matrix times vector (k)->(k,1)
    int m, n, k, ac, bc, cc;
    auto AT = CblasTrans;
    if (aI==0) { // A is "transposed"
      m = A.shape.dims[1];
      k = A.shape.dims[0];
    }
    else {
      m = A.shape.dims[0];
      k = A.shape.dims[1];
      AT = CblasNoTrans;
    }
    ac = A.shape.dims[1];
    auto BT = CblasNoTrans;
    n = bc = cc = 1;
    double ALPHA = 1.0, BETA = 0;

    cblas_dgemm(CblasRowMajor, AT, BT, m, n, k, ALPHA, A.array, ac, B.array, bc, BETA, C.array, cc);
    return;
  }
  
  // STUB
  
  
  
}

void multiply(const Tensor& A, const Tensor& B, Tensor& C) {
  multiply(A, A.shape.rank-1, B, 0, C);
}

void multiply(const double m, const Tensor& A, const Tensor& B) {
  for (int i=0; i<A.total; i++) B.array[i] = m*A.array[i];
}

void timesEq(Tensor& A, const double m) {
  for (int i=0; i<A.total; i++) A.array[i] *= m;
}

void add(const Tensor &A, const Tensor& B, Tensor& C) {
  Tensor::checkDims(A, B); Tensor::checkDims(A, C);
  for (int i=0; i<A.total; i++) C.array[i] = A.array[i] + B.array[i];
}

void NTplusEqUnsafe(Tensor& A, const Tensor& B, double mult) {
  for (int i=0; i<A.total; i++) A.array[i] += mult*B.array[i];
}

void subtract(const Tensor &A, const Tensor& B, Tensor& C) {
  Tensor::checkDims(A, B); Tensor::checkDims(A, C);
  for (int i=0; i<A.total; i++) C.array[i] = A.array[i] - B.array[i];
}

void NTminusEqUnsafe(Tensor& A, const Tensor& B, double mult) {
  for (int i=0;i<A.total; i++) A.array[i] -= mult*B.array[i];
}

void TminusEq(Tensor& A, const Tensor& B, double mult) {
  // Do checks
  if (A.shape.rank!=2 || B.shape.rank!=2) throw Tensor::TensorBadFunction();
  if (A.shape.dims[0]!=B.shape.dims[1] || A.shape.dims[1]!=B.shape.dims[0])
    throw Tensor::TensorDimsMismatch();
  // Do subtraction
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      A.at(i,j) -= mult*B.at(j,i);
}

void hadamard(const Tensor &A, const Tensor& B, Tensor& C) {
  Tensor::checkDims(A, B); Tensor::checkDims(A, C);
  for (int i=0; i<A.total; i++) C.array[i] = A.array[i] * B.array[i];
}

void hadamardEq(Tensor& A, const Tensor& B) {
  Tensor::checkDims(A, B);
  for (int i=0; i<A.total; i++) A.array[i] *= B.array[i];
}

void apply(const Tensor& A, function F, Tensor& C) {
  Tensor::checkDims(A, C);
  for (int i=0; i<A.total; i++) C.array[i] = F(A.array[i]);
}

int Tensor::getDim(int i) {
  if (i<0 || i>shape.rank) throw TensorRankMismatch();
  return shape.dims[i];
}

double Tensor::getSum() {
  double sum = 0;
  for (int i=0; i<total; i++) sum += array[i];
  return sum;
}

void Tensor::resize(const Shape& s) {
  if (array) delete [] array;
  if (stride) delete [] stride;
  *this = Tensor(s);
}

void Tensor::reshape(const Shape& s) {
  int tot = s.getTotal();
  if (total!=tot) throw TensorBadReshape();
  initialize(s, false);
}

void Tensor::random(double max) {
  for (int i=0; i<total; i++)
    array[i] = max*(2*drand48()-1);
}

void Tensor::zero() {
  for (int i=0; i<total; i++) array[i] = 0;
}

void Tensor::qrel() {
  array = 0;
}

void Tensor::qref(Tensor& T) {
  if (array) delete [] array;
  array = T.array;
}

inline void Tensor::writeHelper(vector<int> indices, std::ostream& out, const Tensor& T) {
  out << '{';
  int step = indices.size();
  if (step==T.shape.rank-1) { // Base case
    if (T.shape.dims==0) {
      out << "}";
      return;
    }
    for (int i=0; i<T.shape.dims[step]; i++) {
      indices.push_back(i); // Add i
      out << T.at(indices);
      indices.pop_back();   // Remove i
      if (i!=T.shape.dims[step]-1) out << ',';
    }
  }
  else {
    if (T.shape.dims==0) {
      out << "}";
      return;
    }
    for (int i=0; i<T.shape.dims[step]; i++) {
      indices.push_back(i); // Add i
      writeHelper(indices, out, T);
      indices.pop_back();   // Remove i
      if (i!=T.shape.dims[step]-1) out << ',';
    }
  }
  out << '}';
}

// The stream 'in' should point at the start of the entries, either numbers or subtensors
inline Tensor Tensor::readHelper(std::istream& in) {
  char c; 
  in.get(c);
  if (c=='{') {
   vector<Tensor> entries;   
    while (c!='}') {
      if (c==',') in.get(c);
      Tensor ten = readHelper(in); // Get a rank n subtensor
      entries.push_back(ten);
      if (!in.eof()) in.get(c); // Either get , or }
      else break;
    }
    // Package into a single tensor (rank n+1)
    Tensor T(entries);
    return T;
  }
  else if (isdigit(c)) { // This is the bottom level
    in.putback(c);
    vector<double> entries;
    double n;
    while (c!='}') {
      in >> n;
      entries.push_back(n);
      if (!in.eof()) in.get(c); // Either get , or }
      else break;
    }
    // Package into a single tensor (rank 1, a vector)
    Tensor T(entries);
    return T;
  }
  else return Tensor(); // This shouldn't happen
}

std::ostream& operator<<(std::ostream& out, const Tensor& T) {
  vector<int> indices;
  Tensor::writeHelper(indices, out, T);
  return out;
}

std::istream& operator>>(std::istream& in, Tensor& T) {
  char c;
  vector<int> indices;
  in.get(c);
  // Check for the start of a tensor
  if (c=='{') T = Tensor::readHelper(in);
  else T.reshape(Shape(1)); // Else - Not a tensor
  return in;
}

void Tensor::initialize(Shape s, bool del, bool zero, int tot) {
  shape = s;
  // Find total
  total = s.getTotal();
  // Set stride array
  int count = 1, rank = shape.rank;
  stride = new int[rank];
  for (int i=0; i<rank; i++) {
    count *= shape.dims[i];
    stride[i] = total/count;
  }

  if (del) {
    // Set data array
    array = new double[total];
    if (zero) for (int i=0; i<total; i++) array[i] = 0.;
  }
}

inline bool Tensor::checkDims(const Tensor& A, const Tensor& B) {
  if (A.shape.rank != B.shape.rank) throw TensorRankMismatch();
  for (int i=0; i<A.shape.rank; i++) 
    if (A.shape.dims[i]!=B.shape.dims[i])
      throw TensorDimsMismatch();
}

inline void Tensor::shift_helper(const Shape& shift, vector<int>& point, Tensor& T) const {
  int iter = point.size();
  if (shape.rank-1==iter) {
    for (int i=0; i<shape.dims[iter]; i++) {
      point.push_back(i);
      double value = at(point);
      T.set(value, point, shift);
      point.pop_back();
    }
  }
  else
    for (int i=0; i<shape.dims[iter]; i++) {
      point.push_back(i);
      shift_helper(shift, point, T);
      point.pop_back();
    }
}
