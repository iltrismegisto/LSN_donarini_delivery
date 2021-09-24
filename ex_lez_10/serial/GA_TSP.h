#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

int wd=18;
int num_city=32;
double temp_init;
int temperature_steps;
int steps_per_T;
double temp;
int mode = 0;
double A=0;
int accepted_1 = 0,accepted_2 = 0;

double inline Random(){
    return rand()/(double)RAND_MAX;
}

class City{
  public:
  City();

  City(double x, double y, int i) {
    m_x = x;
    m_y = y;
    index = i;
  }
  ~City(){};

  double GetX() {return m_x;}
  double GetY() {return m_y;}
  int GetIndex() {return index;}

  void SetX(double x) {m_x = x;}
  void SetY(double y) {m_y = y;}
  void SetIndex(int i) {index = i;}

  private:
  double m_x;
  double m_y;
  int index;
};

vector <City> cities;

class Gene{
  public:
  Gene() {
  cout << "--- Empty gene ---" << endl;
  }

  Gene(int length, vector<City>* c, bool shuffle) {
    for(size_t i=0 ; i<length ; i++ ) {
      gene.push_back(c->at(i));
    }
    if(shuffle) Shuffle();
    cost = Cost();
  }
  ~Gene(){
    gene.clear();
  }
  void AddElement(City city) {gene.push_back(city); }

  City GetElement(int i) {return gene[i];}

  int GetLength() {return gene.size();}

  double GetCost() {return Cost(); }

  void SwapCity(int x, int y){
    City temp = City(0,0,0);
    temp = gene[x];
    gene[x] = gene[y];
    gene[y] = temp;
  }

  double Cost() {
    double r = 0.;
    for (size_t i = 0; i < GetLength(); i++) {
      double j = i+1;
      if(i==GetLength()-1) j=0; //PBC
        r += sqrt(pow((gene[i].GetX()-gene[j].GetX()),2)+pow((gene[i].GetY()-gene[j].GetY()),2));
    }
    return r;
  }

  void Shuffle(){
    for (size_t i=0; i<GetLength(); i++){
      random_shuffle(gene.begin()+1, gene.end());
    }
  }

  void Print() {
    cout << "--------cities--------"<<endl;
    for (size_t i=0; i<GetLength(); i++){
      cout << "City #" << gene[i].GetIndex() << " at ("<<  gene[i].GetX() << " , " << gene[i].GetY() << ")"<< endl;
    }
    cout << "Cost: " << Cost() << " $"<< endl;
  }

  void Mutation_Pair(){
    int where = 1+Random()*(GetLength()-2);
    SwapCity(where,where+1);
  }

  void Mutation_Inversion(){
    int how_many = 3;
    int where = 1+Random()*(GetLength()-how_many-1);
    for ( size_t i=0 ; i <= int(how_many/2) ; i++ ){
      SwapCity(where+i, where+how_many-i);
    }
  }

  private:
  int length;
  vector<City> gene;
  double cost;
};

class Population {

  public:
  Population();

  Population(int pop_dim, int num_city) {
    dim=pop_dim;
    for ( size_t i = 0 ; i < dim ; i++ ) {
      Gene gene = Gene(num_city, 0, true);
      population.push_back(gene);
    }
    Sort();
  }

  Population(int pop_dim, int num_city, vector<City>* c){
    dim=pop_dim;
    for (size_t i=0 ; i < dim ; i++) {
      Gene gene = Gene(num_city, c, true);
      population.push_back(gene);
    }
    Sort();
  }

  ~Population() {population.clear();}

  void Clear() {population.clear();}

  void AddGene(Gene gene) {population.push_back(gene);}

  Gene GetGene(int i) {return population[i];}

  int GetDim() {return population.size();}

  static bool sort_by_cost (Gene& x, Gene& y) {return (x.GetCost() < y.GetCost());}

  void Sort() {sort (population.begin(), population.end(), sort_by_cost);}

  void Mutation_Pair(Gene* g){
    int where = 1+Random()*(g->GetLength()-2);
    g->SwapCity(where,where+1);
  }

  void Mutation_Shift(Gene* g){
    int how_many = 3;
    int how_far = how_many+Random()*10;
    int where = 1+Random()*(g->GetLength()-how_far-how_many-1);
    for(int i=0;i<how_many;i++){
      g->SwapCity(where+i, where+how_far+i);
    }
  }

  void Mutation_Permutation(Gene* g){
    int how_many = 3;
    int where_f = 1+Random()*(g->GetLength()-how_many-1);
    int where_s = 1+Random()*(g->GetLength()-how_many-1);
    int where_t = 1+Random()*(g->GetLength()-how_many-1);
    while ( (where_f-how_many) > where_s or where_s > (where_f+how_many) ) {
      where_s = 1+Random()*(g->GetLength()-how_many-1);
    }
    while ( ((where_f-how_many) > where_t or where_t > (where_f+how_many)) and ((where_s-how_many) > where_t or where_t > (where_s+how_many))) {
      where_t = 1+Random()*(g->GetLength()-how_many-1);
    }

    for (int i=0 ; i < how_many ; i++ ) {
      g->SwapCity(where_f+i, where_s+i );
      g->SwapCity(where_s+i, where_t+i );
    }
  }

  void Mutation_Inversion(Gene* g){
    int how_many = 3;
    int where = 1+Random()*(g->GetLength()-how_many-1);
    for ( size_t i=0 ; i <= int(how_many/2) ; i++ ){
      g->SwapCity(where+i, where+how_many-i);
    }
  }

  void Mutate(double mutation_rate){
    for(size_t i=0; i<GetDim(); i++){
      if(Random()<mutation_rate) {
        //cout << "Pair"<< endl;
        Mutation_Pair(&population[i]);

      }
      if(Random()<mutation_rate){
        //cout << "Shift"<< endl;
        Mutation_Shift(&population[i]);

      }
      if(Random()<mutation_rate){
        //cout << "Permutation"<< endl;
        Mutation_Permutation(&population[i]);

      }
      if(Random()<mutation_rate){
        //cout << "Inversion"<< endl;
        Mutation_Inversion(&population[i]);

      }
    }
  }

  Gene SelectGene(double exp){
    double s = Random();
    double prob = pow(s,exp);
    double pos = int(dim*prob); // manca il +1 perché non capisco mi dia sempre segmentation fault
    return population[pos];
  }

  void Cross(Gene x, Gene y, Population* new_pop) {
  int length = x.GetLength();
  int where = 1+Random()*length;
  Gene a = Gene(0,0,false);
  Gene b = Gene(0,0,false);
  for ( size_t i=0 ; i < where ; i++ ) {
    a.AddElement(x.GetElement(i));
    b.AddElement(y.GetElement(i));
  }
  while( a.GetLength() < length){
    for( size_t i=0 ; i < length ; i++ ){
      bool nw = true;
      for (size_t j=0; j<a.GetLength(); j++) {
        if(y.GetElement(i).GetIndex() == a.GetElement(j).GetIndex()) nw = false;
      }
      if ( nw ) {
        a.AddElement(y.GetElement(i));
      }
    }
  }
  while(b.GetLength()<length){
    for( size_t i=0 ; i<length ; i++ ){
      bool nw = true;
      for (size_t j=0; j < b.GetLength(); j++)
      {
        if(x.GetElement(i).GetIndex()==b.GetElement(j).GetIndex()) nw=false;
      }
      if( nw ){
        b.AddElement(x.GetElement(i));
      }
    }
  }
  new_pop->AddGene(a);
  new_pop->AddGene(b);
  }

  void SaveBest(int saved_per_generation, double selection_exp, Population* new_pop){
    for (int i = 0; i < saved_per_generation; i++){
    Gene a = SelectGene(selection_exp);
    new_pop->AddGene(a);
    }
  }

  void Print(){
    cout << "--------Population--------"<<endl;
    for (size_t i=0; i<GetDim(); i++){
      cout << "Gene #" << i << " | Cost: "  <<population[i].GetCost() << endl;
      }
    }
  private:
  int dim;
  vector<Gene> population;
};

void Input( vector<City> *cities ) {

  cout << "---START SIMULATED ANNEALING---" << endl;
  ifstream input;
  input.open("input.dat");

  input >> num_city;
  cout << "Num city : " << num_city << endl;

  input >> temp_init;
  cout << "Initial temperature : " << temp_init << endl;
  input >> temperature_steps;
  cout << "Temperature steps = " << temperature_steps << endl;
  input >> steps_per_T;
  cout << "Steps at fixed temperature = " << steps_per_T << endl;
  input >> mode;
  if (mode==0) {
    cout << "Mode = Circle" << endl << endl;
  } else {
    cout << "Mode = Square" << endl << endl;
  }
  input.close();

  if(mode == 0) {   //Cities on a circumference
    double angle = 2*M_PI/double(num_city);
    cout << "Cities generated on a circumference "<< endl << endl;
    for (size_t i=0; i<num_city; i++){
      City city = City(cos(angle*i), sin(angle*i) ,i);
      cities->push_back(city);
    }
  }

  if(mode == 1) {   //Cities inside a square
    cout << "Cities generated inside a square "<< endl << endl;
    for (size_t i=0; i<num_city; i++){
      City city = City(Random(), Random() ,i);
      cities->push_back(city);
    }
  }

}
