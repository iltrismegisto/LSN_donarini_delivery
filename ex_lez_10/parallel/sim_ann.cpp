#include "mpi.h"
#include "GA_TSP_PAR.h"
#include <sstream>

using namespace std;

int main(int argc, char* argv[]){
  MPI::Init(argc,argv);
  int size = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  srand(rank);

  Input(rank);

  Population pop = Population(dim, num_city, &cities);
  Population next_pop = Population(0, num_city, &cities);

  ofstream output, best, first, out_best_rank;

  double node_cost_array[size];
  double node_cost=0;
  int best_rank = 0;
  double best_gene_cost = 0;
  int n_mig = Random()*num_generations;

  for(size_t i=0;i<size;i++) node_cost_array[i]=0;

  Gene g = pop.GetGene(0);

  if (rank ==0 ) {
    output.open("evolution_sq.dat");
    out_best_rank.open("best_rank_evo.dat");

  }


  for (size_t ranks = 0; ranks < size; ranks++ ) {

    if (rank == ranks) {
      std::ostringstream fn;
      fn << "first_sq_" << ranks << ".dat";
      first.open(fn.str().c_str());
      

      for (size_t i = 0; i < g.GetLength(); i++){
        City c = g.GetElement(i);
        first << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;

      }
    first.close();
    }
  }
  
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

    if (rank==0 )cout << "Generation #" << i << " best Salesman Cost: " << pop.GetGene(0).GetCost() <<endl;
    if (rank==0 )output << setw(wd) << i <<  setw(wd) << pop.GetGene(0).GetCost() << endl;
    
    next_pop.Clear();

    node_cost = pop.GetGene(0).GetCost();

    MPI_Gather(&node_cost,1,MPI_DOUBLE_PRECISION,node_cost_array,1,MPI_DOUBLE_PRECISION,0,MPI::COMM_WORLD);

    
    if(rank==0){
      best_gene_cost = 999999;
        for (int k = 0; k < size; k++){
            if(node_cost_array[k]<best_gene_cost){
                best_gene_cost = node_cost_array[k];
                best_rank = k;
            }
        }
    }
    if (rank==0) out_best_rank << setw(wd) << i <<  setw(wd) << best_rank << endl; 

    if (num_generations==n_mig){
        int w1 = Random()*4;
        int w2 = Random()*4;
        bool wait=false;
        int howmany = saved_per_generation/2;
        Population pop_temp=Population(0, num_city, &cities);
        if (w2==w1) {
           w2 = Random()*4;
        }

    	if (rank==w1){

            for (size_t y=0; y<dim; y++){
                pop_temp.AddGene(pop.GetGene(y));
            }
            MPI_Recv(&pop, 1, MPI_BYTE, w2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            wait=true;
            }
            if (rank == w2) {
                while (wait==true) {
                }
                for (size_t y=0; y<howmany; y++){
                    pop.AddGene(pop_temp.GetGene(y));
                    pop.Sort();
                }
            }
        pop_temp.Clear();
      }
    next_pop.Clear();
  }

  output.close();

  MPI_Bcast(&best_rank, 1, MPI_INTEGER, 0, MPI::COMM_WORLD);

  if (rank==best_rank){

      best.open("best_sq.gene");
      g = pop.GetGene(0);
      for (size_t i = 0; i < g.GetLength(); i++){
        City c = g.GetElement(i);
        best << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
      }

  best.close();
  MPI::Finalize();

  return 0;
}

  MPI::Finalize();
  return 0;

}

