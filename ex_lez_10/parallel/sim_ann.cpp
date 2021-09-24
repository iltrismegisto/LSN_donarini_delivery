#include "mpi.h"
#include "GA_TSP.h"

using namespace std;

int main(int argc, char* argv[]){
  MPI::Init(argc,argv);
  int size = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  srand(rank);

  Input(rank, &cities);

  ofstream output, best, first, out_best_rank;

  Gene gene = Gene(num_city,&cities,true);
  Gene new_gene = Gene(num_city,&cities,true);
  Gene best_gene = gene;

  double node_cost_array[size];
  double node_cost=0;
  int best_rank = 0;
  double best_gene_cost = 0;

  for(size_t i=0;i<size;i++) node_cost_array[i]=0;

  if(rank==0){

    if (temp_init==1) {
      if(mode == 0)output.open("evolution_t1.dat");
      if(mode == 1)output.open("evolution_t1_sq.dat");
      if(mode == 0)first.open("first_t1.gene");
      if(mode == 1)first.open("first_t1_sq.gene");
      if(mode == 0)out_best_rank.open("best_rank_evo_t1.dat");
      if(mode == 1)out_best_rank.open("best_rank_evo_t1_sq.dat");
    }

    if (temp_init==3) {
      if(mode == 0)output.open("evolution_t2.dat");
      if(mode == 1)output.open("evolution_t2_sq.dat");
      if(mode == 0)first.open("first_t2.gene");
      if(mode == 1)first.open("first_t2_sq.gene");
      if(mode == 0)out_best_rank.open("best_rank_evo_t2.dat");
      if(mode == 1)out_best_rank.open("best_rank_evo_t2_sq.dat");
    }

    if (temp_init==10) {
      if(mode == 0)output.open("evolution_t3.dat");
      if(mode == 1)output.open("evolution_t3_sq.dat");
      if(mode == 0)first.open("first_t3.gene");
      if(mode == 1)first.open("first_t3_sq.gene");
      if(mode == 0)out_best_rank.open("best_rank_evo_t3.dat");
      if(mode == 1)out_best_rank.open("best_rank_evo_t3_sq.dat");
    }

    for (size_t i = 0; i < gene.GetLength(); i++){
      City c = gene.GetElement(i);
      first << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
    }
    first.close();
  }

  for (size_t i = temperature_steps; i > 0; i--){
    temp = (temp_init/(double)temperature_steps)*i;
    if(i%50==0 && rank==0)cout<< "---------- Temperature " << temp << " ----------" << endl;
    accepted_1 = 0;
    accepted_2 = 0;
    for (size_t j = 0; j < steps_per_T; j++) {
      new_gene = gene;
      new_gene.Mutation_Pair();
      A = min(1.,exp(-(new_gene.Cost()-gene.Cost())/temp));
      if(Random()<A) {
        gene = new_gene;
        accepted_1 ++;
      }
      if(best_gene.Cost()>gene.Cost()) best_gene=gene;

      new_gene = gene;
      new_gene.Mutation_Inversion();
      A = min(1.,exp(-(new_gene.Cost()-gene.Cost())/temp));
      if(Random()<A) {
        gene = new_gene;
        accepted_2 ++;}
      if(best_gene.Cost()>gene.Cost()) best_gene=gene;

    }
    node_cost = best_gene.Cost();
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

    if(rank==0) output << setw(wd) << temperature_steps-i <<  setw(wd) << best_gene_cost << endl;
    if(rank==0) out_best_rank << setw(wd) << temperature_steps-i <<  setw(wd) << best_rank << endl;
    if(i%50==0 && rank==0) {
      cout<< "Pair Acceptation rate at temp " << temp << " is: " << (double)accepted_1/(double)steps_per_T << endl;
      cout<< "Inversion Acceptation rate at temp " << temp << " is: " << (double)accepted_2/(double)steps_per_T << endl;
      cout<< "Best gene at temp " << temp << " has cost: " << best_gene_cost << endl << endl;
    }
  }



  MPI_Bcast(&best_rank, 1, MPI_INTEGER, 0, MPI::COMM_WORLD);

  if(rank==0) cout<< "**** Best gene has cost: " << best_gene_cost << endl << endl;


  if (rank==best_rank) {
    if (temp_init==1) {
      if(mode == 0)best.open("best_t1.gene");
      if(mode == 1)best.open("best_t1_sq.gene");
    }

    if (temp_init==3) {
      if(mode == 0)best.open("best_t2.gene");
      if(mode == 1)best.open("best_t2_sq.gene");
    }

    if (temp_init==10) {
      if(mode == 0)best.open("best_t3.gene");
      if(mode == 1)best.open("best_t3_sq.gene");
    }

    for (size_t i = 0; i < best_gene.GetLength(); i++){
        City c = best_gene.GetElement(i);
        best << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
    }
  }
  best.close();




  if(rank==0) {
    output.close();
    out_best_rank.close();
  }

  MPI::Finalize();
  return 0;

}
