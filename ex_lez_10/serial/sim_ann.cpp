#include "GA_TSP.h"

using namespace std;

int main(){
    Input(&cities);

    ofstream output, best, first;

    if (temp_init==1) {
      if(mode == 0)output.open("evolution_t1.dat");
      if(mode == 1)output.open("evolution_t1_sq.dat");
      if(mode == 0)first.open("first_t1.gene");
      if(mode == 1)first.open("first_t1_sq.gene");
      if(mode == 0)best.open("best_t1.gene");
      if(mode == 1)best.open("best_t1_sq.gene");
    }

    if (temp_init==3) {
      if(mode == 0)output.open("evolution_t2.dat");
      if(mode == 1)output.open("evolution_t2_sq.dat");
      if(mode == 0)first.open("first_t2.gene");
      if(mode == 1)first.open("first_t2_sq.gene");
      if(mode == 0)best.open("best_t2.gene");
      if(mode == 1)best.open("best_t2_sq.gene");
    }

    if (temp_init==10) {
      if(mode == 0)output.open("evolution_t3.dat");
      if(mode == 1)output.open("evolution_t3_sq.dat");
      if(mode == 0)first.open("first_t3.gene");
      if(mode == 1)first.open("first_t3_sq.gene");
      if(mode == 0)best.open("best_t3.gene");
      if(mode == 1)best.open("best_t3_sq.gene");
    }

    Gene gene = Gene(num_city,&cities,true);
    Gene new_gene = Gene(num_city,&cities,true);
    Gene best_gene = gene;

    for (size_t i = 0; i < gene.GetLength(); i++){
        City c = gene.GetElement(i);
        first << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
    }
    first.close();

    for (size_t i = temperature_steps; i > 0; i--){
        temp = (temp_init/double(temperature_steps))*i;
        if(i%50==0) cout<< "----- Temperature " << temp << " ------" << endl;
        accepted_1 = 0;accepted_2 = 0;
        for (size_t j = 0; j < steps_per_T; j++) {
            new_gene = gene;
            new_gene.Mutation_Pair();
            A = min(1.,exp(-(new_gene.Cost()-gene.Cost())/temp));
            if(Random()<A) { gene = new_gene; accepted_1 ++;}
            if(best_gene.Cost()>gene.Cost()) best_gene=gene;

            new_gene = gene;
            new_gene.Mutation_Inversion();
            A = min(1.,exp(-(new_gene.Cost()-gene.Cost())/temp));
            if(Random()<A){ gene = new_gene; accepted_2 ++;}
            if(best_gene.Cost()>gene.Cost())best_gene=gene;

        }
        output << setw(wd) << temperature_steps-i <<  setw(wd) << best_gene.GetCost() << endl;
        if(i%50==0){
            cout<< "Pair Acceptation rate at temp " << temp << " is: " << (double)accepted_1/(double)steps_per_T << endl;
            cout<< "Inversion Acceptation rate at temp " << temp << " is: " << (double)accepted_2/(double)steps_per_T << endl;
            cout<< "Best gene at temp " << temp << " has cost: " << best_gene.Cost() << endl << endl;
        }
    }

    cout<< "Best gene: " << best_gene.Cost() << endl << endl;
    for (size_t i = 0; i < best_gene.GetLength(); i++){
        City c = best_gene.GetElement(i);
        best << setw(wd) << c.GetIndex() <<  setw(wd) << c.GetX() << setw(wd) << c.GetY() << endl;
    }

    best.close();

    return 23;
}
