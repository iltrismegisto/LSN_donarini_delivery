#include "GA_TSP.h"

using namespace std;

int main(){

  Input();

  Population pop = Population(dim, num_city, &cities);
  Population next_pop = Population(0, num_city, &cities);

  ofstream output, best, first;
  if(mode == 0)output.open("evolution.dat");
  if(mode == 1)output.open("evolution_sq.dat");
  if(mode == 0)first.open("first.gene");
  if(mode == 1)first.open("first_sq.gene");

  Gene g = pop.GetGene(0);
  for (size_t i = 0; i < g.GetLength(); i++){
    City c = g.GetElement(i);
    first << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;

  }
  first.close();

  for (size_t i = 0; i < num_generations; i++) {

    pop.SaveBest(saved_per_generation, selection_exp, &next_pop);
    pop.Mutate(mutation_rate);
    pop.Sort();

    if((dim-saved_per_generation)%2!=0) {
      Gene x = pop.SelectGene(selection_exp);
      next_pop.AddGene(x);
    }










    
    for(size_t j = 0; j < int(dim-saved_per_generation)/2; j++){
      Gene x = pop.SelectGene(selection_exp);
      Gene y = pop.SelectGene(selection_exp);
      if(Random()<cross_rate) {
        pop.Cross(x,y, &next_pop);
      }
      else{
        next_pop.AddGene(x);next_pop.AddGene(y);
      }
    }











    next_pop.Sort();
    pop.Clear();
    for (size_t j = 0; j < next_pop.GetDim(); j++){
      pop.AddGene(next_pop.GetGene(j));
    }

    pop.Mutate(mutation_rate);
    pop.Sort();

    cout << "Generation #" << i << " best Salesman Cost: " << pop.GetGene(0).GetCost() <<endl;
    output << setw(wd) << i <<  setw(wd) << pop.GetGene(0).GetCost() << endl;
    next_pop.Clear();


  }


  output.close();
  if(mode == 0) best.open("best.gene");
  if(mode == 1) best.open("best_sq.gene");
  g = pop.GetGene(0);
  for (uint i = 0; i < g.GetLength(); i++){
    City c = g.GetElement(i);
    best << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
  }
  best.close();

  return 23;
}
