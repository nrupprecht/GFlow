#include "GField.h"

GField::GField() : needsRemake(false){
  zeroPointers();
}

GField::GField(int rk) {
  rank = rk;
  zeroPointers();
  init(rk, true);
}

GField::GField(Shape s) : needsRemake(false) {
  zeroPointers();
  initialize(s, true);
}

GField::GField(int* sizes, int dims) : needsRemake(false){
  zeroPointers();
  Shape s(sizes, dims);
  initialize(s, true);
}

GField::~GField() {
  clean();
}

GField& GField::operator=(const GField& F) {
  clean();
  shape = F.shape;
  total = F.total;
  rank = F.rank;
  needsRemake = F.needsRemake;
  create(F.spacing, spacing, rank);
  create(F.stride, stride, rank);
  create(F.array, array, total);
  create(F.lowerBounds, lowerBounds, rank);
  create(F.upperBounds, upperBounds, rank);
  create(F.wrapping, wrapping, rank);
}

inline void GField::at_address(int& add, Index index) const { // Doesn't check (for speed)
  for (int i=0; i<index.size(); i++) add += stride[i]*index.at(i);
}

inline void GField::at_address(int& add, vector<int> index) const { // Doesn't check (for speed)
  for (int i=0; i<index.size(); i++) add += stride[i]*index.at(i);
}

void GField::initialize(Shape s, bool all) {
  int lastRank = shape.getRank(), lastTotal = total;
  if (s.getRank()==0 || s.getTotal()==0) return;
  total = s.getTotal();
  shape = s;
  // Initialize other arrays
  init(s.getRank(), all);
  // Set up array (recalculate on change of Shape)
  if (total!=lastTotal || all) {
    if (array) delete [] array;
    array = nullptr;
  }
  if (array==nullptr) array = new double[total];
  for (int i=0; i<total; i++) array[i] = 0.; // Always reset values
}

double& GField::at(Index index) {
  int address = 0;
  at_address(address, index);
  if (total<=address) throw GFieldOutOfBounds();
  return array[address];
}

double GField::at(Index index) const {
  int address = 0;
  at_address(address, index);
  if (total<=address) throw GFieldOutOfBounds();
  return array[address];
}

double& GField::at(vector<int> index) {
  int address = 0;
  at_address(address, index);
  if (total<=address) throw GFieldOutOfBounds();
  return array[address];
}

double GField::at(vector<int> index) const {
  int address = 0;
  at_address(address, index);
  if (total<=address) throw GFieldOutOfBounds();
  return array[address];
}

double GField::derivative(Index I, Index pos) const {
  // Assumes a single first derivative
  Index up = pos+I, lw = pos-I;
  int ind1 = mod(up);
  int ind2 = mod(lw);
  // Find which index we are differentiating
  int index = 0;
  for (int i=0; i<I.size(); i++)
    if (I.at(i)==1) {
      index = i;
      break;
    }
  double invh = 1./spacing[index];
  double der = 0;
  if (ind1==0 && ind2==0) der = (at(up)-at(lw))*0.5*invh;
  else if (ind1!=0) der = (at(I)-at(lw))*invh;
  else der = (at(up)-at(I))*invh;

  return der;
}

gVector GField::getPosition(vector<int> indices) const {
  gVector pos;
  pos.resize(indices.size());
  for (int i=0; i<indices.size(); i++) 
    pos.at(i) = indices.at(i)*spacing[i]+lowerBounds[i];
  return pos;
}

gVector GField::getPosition(Index indices) const {
  gVector pos;
  pos.resize(indices.size());
  for (int i=0; i<indices.size(); i++)
    pos.at(i) = indices.at(i)*spacing[i]+lowerBounds[i];
  return pos;
}

gVector GField::getPosition(int i) const {
  Index index = getIndex(i); // Find the index
  return getPosition(index); // Use the index get position function
}

double GField::getPosition(int index, int grid) const {
  if (shape.getRank()<=index) throw GFieldIndexOutOfBounds();
  return grid*spacing[index]+lowerBounds[index];
}

Index GField::getIndex(int i) const {
  int num = total;
  vector<int> indices;
  for (int j=0; j<shape.getRank(); j++) {
    num /= shape.at(j);
    int I = j/num;
    indices.push_back(I);
    i %= num;
  }
  return Index(indices);
}

double GField::integrate() {
  double total = 0.;
  for (int i=0; i<shape.getTotal(); i++) total += array[i];
  double dV = 1;
  for (int i=0; i<shape.getRank(); i++) dV *= spacing[i];
  return total*dV;
}

double GField::integrate(int index, Index pos) {
  int add = 0;
  pos.at(index) = 0; // Start at the beginning
  at_address(add, pos);
  // Integrate
  int step = stride[index];
  double total = 0.;
  for (int i=0; i<shape.at(index); i++, add+=step) total += array[add];  
  total *= spacing[index];
  return total;
}

GField& GField::operator+=(const GField& F) {
  if (F.shape!=shape) throw GFieldMismatch();
  for (int i=0; i<total; i++)
    array[i] += F.array[i];
}

/// Not safe (no checking) version
void plusEqNS(GField& G, const GField& H, double mult) {
  for (int i=0; i<G.total; i++) G.array[i] += mult*H.array[i];
}

int GField::getSize(int index) const {
  if (index<0 || shape.getRank()<=index) throw GFieldIndexOutOfBounds();
  return shape.at(index);
}

int GField::getPoints() const {
  if (shape.getRank()==0) return 0;
  int points = 1;
  for (int i=0; i<shape.getRank(); i++) points *= shape.at(i);
  return points;
}

void GField::setWrapping(int i, bool w) {
  if (i<0 || rank<=i) GFieldDimensionMismatch();  
  wrapping[i] = w;
}

void GField::setBounds(int i, double l, double u, bool remake) {
  if (i<0 || rank<=i) throw GFieldDimensionMismatch();
  upperBounds[i] = u; lowerBounds[i] = l;
  // Recalculate the spacing for this dimension
  if (shape!=Shape())
    spacing[i] = (upperBounds[i]-lowerBounds[i])/shape.at(i);
  // Remake not needed, clear array if requested
  if (remake) for (int i=0; i<shape.getTotal(); i++) array[i] = 0;
  else needsRemake = true;
}

void GField::setGridSpacing(int i, double delta, bool remake) {
  if (i<0 || rank<=i) throw GFieldDimensionMismatch();
  spacing[i] = delta;
  int size = (upperBounds[i]-lowerBounds[i])/spacing[i];
  if (shape!=Shape()) shape.set(i, size);
  // Remake if requested
  if (remake) initialize(shape, false);
  else needsRemake = true;
}

void GField::set(gFunction function) {
  Index index;
  index.resize(shape.getRank());
  GField::setHelper(index, 0, function, *this);
}

void GField::parabolaZero(int index, int N) {
  //** STUB
}

void GField::remake() {
  if (shape==Shape()) createShape();
  initialize(shape, false);
  needsRemake = false;
}

std::ostream& operator<<(std::ostream& out, const GField& G) {
  vector<int> indices;
  string str = "{", str2;
  stringstream stream;
  GField::writeHelper(indices, stream, G);
  stream >> str2;
  str += str2;
  str.pop_back(); // Put back the last ','
  str += "}";
  out << str;
  return out;
}

inline void GField::writeHelper(vector<int> indices, std::ostream& out, const GField& G) {
  int step = indices.size();
  if (step==G.shape.getRank()-1) { // Base case
    for (int i=0; i<G.shape.at(step); i++) {
      indices.push_back(i); // Add i
      gVector coord = G.getPosition(indices);
      out << '{' << bare_representation(coord) << ',' << G.at(indices) << '}';
      indices.pop_back();   // Remove i
      out << ",";
    }
  }
  else {
    for (int i=0; i<G.shape.at(step); i++) {
      indices.push_back(i); // Add i
      writeHelper(indices, out, G);
      indices.pop_back();   // Remove i
    }
  }
}

std::istream& operator>>(std::ostream& in, GField& G) {
  //** STUB
}

inline void GField::setHelper(Index index, int step, gFunction function, GField& G) {
  if (step==G.shape.getRank()-1) { // Base case
    index.at(step)=0;
    for (int i=0; i<G.shape.at(step); i++) {
      gVector coord = G.getPosition(index);
      G.at(index) = function(coord);
      index.at(step)++;
    }
  }
  else {
    index.at(step)=0;
    for (int i=0; i<G.shape.at(step); i++) {
      setHelper(index, step+1, function, G);
      index.at(step)++;
    }
  }
}

inline void GField::init(int rk, bool all) {
  // Set rank
  int lastRank = rank;
  rank = rk;
  // Set stride array (recalculate on change of Shape)
  if (stride) delete [] stride;
  stride = new int[rank];
  if (shape!=Shape()) calculateStride();
  // Set up bounds (recalculate when given)
  if (rank!=lastRank || all) {
    if (lowerBounds) delete [] lowerBounds;
    lowerBounds = new double[rank];
    for (int i=0; i<rank; i++) lowerBounds[i] = -1.;
    if (upperBounds) delete [] upperBounds;
    upperBounds = new double[rank];
    for (int i=0; i<rank; i++) upperBounds[i] = 1.;
  }
  // Set up wrapping (recalculate when given)
  if (rank!=lastRank || all) {
    if (wrapping) delete [] wrapping;
    wrapping = new bool[rank];
    for (int i=0; i<rank; i++) wrapping[i] = false;
  }
  // Set up spacing (recalculate on change of Shape or bounds)
  if (spacing) delete [] spacing;
  spacing = new double[rank];
  if (shape!=Shape()) calculateSpacing();
}

int GField::mod(Index& index) const {
  int oob = 0;
  for (int i=0; i<index.size(); i++) {
    int &c = index.at(i);
    bool w = wrapping[i];
    if (c<0) {
      if (w) c += shape.at(i);
      else oob = -1;
    }
    else if (shape.at(i)<=c) {
      if (w) c -= shape.at(i);
      else oob = 1;
    }
  }
  return oob;
}

inline void GField::clean() {
  if (spacing) delete [] spacing;
  if (stride) delete [] stride;
  if (array) delete [] array;
  if (lowerBounds) delete [] lowerBounds;
  if (upperBounds) delete [] upperBounds;
  if (wrapping) delete [] wrapping;
  zeroPointers();
}

inline void GField::createShape() {
  vector<int> dims;
  for (int i=0; i<rank; i++) 
    dims.push_back((upperBounds[i]-lowerBounds[i])/spacing[i]);
  shape = Shape(dims);
}

inline void GField::zeroPointers() {
  spacing = nullptr;
  stride = nullptr;
  array = nullptr;
  lowerBounds = nullptr;
  upperBounds = nullptr;
  wrapping = nullptr;
}

inline void GField::calculateSpacing() {
  for (int i=0; i<shape.getRank(); i++)
    spacing[i] = (upperBounds[i]-lowerBounds[i])/shape.at(i);
}

inline void GField::calculateStride() {
  int count = 1;
  for (int i=0; i<rank; i++) {
    count *= shape.at(i);
    stride[i] = total/count;
  }
}
