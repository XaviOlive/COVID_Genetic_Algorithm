#include "Fit_documentation.h"

//CONSTANTS
#define popsize 200
#define IC_GENES_NUMBER 3
#define PARAMETERS_GENES_NUMBER 11
#define bits_per_individual 280


void ExitError(const char *miss, int errcode) { //Exit error function
fprintf (stderr, "\nERROR: %s.\nStopping...\n\n", miss); exit(errcode);
}

double uniform(){
    double x;
    return (double)rand() / (double)(RAND_MAX);
} //Function to generate a the pseudorandom numbers between 0 and 1 for mutations and the crossover

unsigned long random_long(unsigned long min, unsigned long max){
    return min + rand() % (max+1 - min);
} //Function to generate the pseudorandom numbers for the initial numbers


void set_individuals(individual * p){
    p->fitness = 0.0;
    p->IC[0] = random_long(0, pow(2,30));
    p -> IC[1] = random_long(0, pow(2,30));
    p -> IC[2] = random_long(0, pow(2,30));
    p -> Pars[0] = random_long(0, pow(2,40));
    p -> Pars[1] = random_long(0, pow(2,20));
    p -> Pars[2] = random_long(0, pow(2,10));
    p -> Pars[3]= random_long(0, pow(2,10));
    p -> Pars[4]= random_long(0, pow(2,20));
    p -> Pars[5]= random_long(0, pow(2,20));
    p -> Pars[6]= random_long(0, pow(2,20));
    p -> Pars[7]= random_long(0, pow(2,10));
    p -> Pars[8]= random_long(0, pow(2,20));
    p -> Pars[9]= random_long(0, pow(2,40));
    p -> Pars[10]= random_long(0, pow(2,40));
}; //Funtion to generate the inituial parametres for an individual



individual* fittest_ind(individual* generation) {
    individual* fittest = &generation[0];
  for (int i = 1; i < popsize; i++) if (generation[i].fitness < fittest->fitness) fittest = &generation[i];
  return fittest;
} //Function that allows us to identify the best individual in the current population


void printIndividual(individual * ind) {
    printf("\n");
  for (int j= 0; j < IC_GENES_NUMBER; j++) {
    printf("IC %d: %lu  -  %f\n", j, ind->IC[j], crom2IC(ind->IC[j]));
  }
    printf("Parametre 1: %lu  -  %.12f\n", ind->Pars[0], crom2HSPar(ind->Pars[0]));
    printf("Parametre 2: %lu  -  %.6f\n", ind->Pars[1], crom2Par(ind->Pars[1]));
    printf("Parametre 3: %lu  -  %.3f\n", ind->Pars[2], crom2LSPar(ind->Pars[2]));
    printf("Parametre 4: %lu  -  %.3f\n", ind->Pars[3], crom2LSPar(ind->Pars[3]));
    printf("Parametre 5: %lu  -  %.6f\n", ind->Pars[4], crom2Par(ind->Pars[4]));
    printf("Parametre 6: %lu  -  %.6f\n", ind->Pars[5], crom2Par(ind->Pars[5]));
    printf("Parametre 7: %lu  -  %.6f\n", ind->Pars[6], crom2Par(ind->Pars[6]));
    printf("Parametre 8: %lu  -  %.3f\n", ind->Pars[7], crom2LSPar(ind->Pars[7]));
    printf("Parametre 9: %lu  -  %.6f\n", ind->Pars[8], crom2Par(ind->Pars[8]));
    printf("Parametre 10: %lu  -  %.12f\n", ind->Pars[9], crom2HSPar(ind->Pars[9]));
    printf("Parametre 11: %lu  -  %.12f\n", ind->Pars[10], crom2HSPar(ind->Pars[10]));
  printf("fitness: %lf\n\n", ind->fitness);
} //Funnction to print the individual (it prints the calculation figure and the real value it represents)


individual* TournamentSelection(individual* population) {
  int decA = uniform() * (popsize - 1), decB = uniform() * (popsize - 1);
  while (decB == decA) decB = uniform() * (popsize - 1);
  if (population[decA].fitness < population[decB].fitness) return (&population[decA]);
  return (&population[decB]);
} //Tournament selection function (it takes the 2 individuals that will be used to generate 2 childs for the new generation)


void OnePointCrossover(individual* predecesor1, individual* predecesor2, individual* child1, individual* child2) {
    int interval_term = (bits_per_individual) * uniform();
    unsigned long provisionalp1[IC_GENES_NUMBER+PARAMETERS_GENES_NUMBER],  provisionalp2[IC_GENES_NUMBER+PARAMETERS_GENES_NUMBER], provisionalc1[IC_GENES_NUMBER+PARAMETERS_GENES_NUMBER], provisionalc2[IC_GENES_NUMBER+PARAMETERS_GENES_NUMBER];
    
    for (int j= 0; j < IC_GENES_NUMBER; j++) {
        provisionalp1[j] = predecesor1->IC[j];
        provisionalp2[j] = predecesor2->IC[j];
    }
    for(int i = 0; i<PARAMETERS_GENES_NUMBER; i++) {
        provisionalp1[i+IC_GENES_NUMBER] = predecesor1->Pars[i];
        provisionalp2[i+IC_GENES_NUMBER] = predecesor2->Pars[i];
    }
    
    unsigned short change = interval_term%14;
    for (int i = 0; i < change; i++) { provisionalc1[i] = provisionalp1[i]; provisionalc2[i] = provisionalp2[i]; }
    for (int i = change; i < IC_GENES_NUMBER+PARAMETERS_GENES_NUMBER; i++) { provisionalc1[i] = provisionalp2[i]; provisionalc2[i] = provisionalp1[i]; }
    
    for(int j= 0; j < IC_GENES_NUMBER; j++) child1->IC[j] = provisionalc1[j];
    for(int i = 0; i<PARAMETERS_GENES_NUMBER; i++) child1->Pars[i] = provisionalc1[i+IC_GENES_NUMBER];
    for(int j= 0; j < IC_GENES_NUMBER; j++) child2->IC[j] = provisionalc2[j];
    for(int i = 0; i<PARAMETERS_GENES_NUMBER; i++) child2->Pars[i] = provisionalc2[i+IC_GENES_NUMBER];
} //ONE_POINT_CROSSOVER FUNCTION: given 2 individuals it pseudorandmly selects a crossover point and 2 chils will we created with the cutted in crossover parametres from one parent and the other.


void Mutation(individual* candidate) {
    int interval_term, interval_term2;
    unsigned long provisional[IC_GENES_NUMBER+PARAMETERS_GENES_NUMBER];
    
    for (int j= 0; j < IC_GENES_NUMBER; j++) provisional[j] = candidate->IC[j];
    for(int i = 0; i<PARAMETERS_GENES_NUMBER; i++) provisional[i+IC_GENES_NUMBER] = candidate->Pars[i];

    interval_term = (bits_per_individual) * uniform();
    interval_term2 = (bits_per_individual) * uniform();
    int j;
    j = interval_term%14;
    if ((j==0) || (j==1) || (j==2)) {
        int i = interval_term%24;
        int k = interval_term2%24;
        provisional[j] = (provisional[j])^(6U << (24-i));
        provisional[j] = (provisional[j])^(6U << (24-k));
    }
    else if ((j==3) || (j==12) || (j==13)) {
        int i = interval_term%30;
        int k = interval_term2%30;
        provisional[j] = (provisional[j])^(10U << (30-i));
        provisional[j] = (provisional[j])^(10U << (30-k));
    }
    else if ((j==5) || (j==6) || (j==10)) {
        int i = interval_term%8;
        int k = interval_term2%8;
        provisional[j] = (provisional[j])^(2U << (8-i));
        provisional[j] = (provisional[j])^(2U << (8-k));
    }
    else {
        int i = interval_term%16;
        int k = interval_term2%16;
        provisional[j] = (provisional[j])^(4U << (16-i));
        provisional[j] = (provisional[j])^(4U << (16-k));
    }
  
    for (int j= 0; j < IC_GENES_NUMBER; j++) candidate->IC[j] = provisional[j];
    for(int i = 0; i<PARAMETERS_GENES_NUMBER; i++)candidate->Pars[i] = provisional[i+IC_GENES_NUMBER];
}// MUTATION FUNCTION: It selects a single parametre of the individual and introduces in its genotype 2 mutations. Dependending on the type of parametres (with different factor from genotype to phenotype) the mutation will afect more or less digits.




int main() {
    // Declarations of variables used and its memory storage (malloc):
    individual * population = NULL;
    DataForFitting * DataPredict = NULL;
    if ((population = (individual*) malloc(popsize*sizeof(individual))) == NULL) ExitError("could not allocate memory for the population", 1);
    if ((DataPredict = (DataForFitting*) malloc(sizeof(DataForFitting))) == NULL) ExitError("could not allocate memory for the DataPredict", 1);
    
    //Definition of the original data given
    DataPredict->PopSize = 1000000;
    DataPredict->N_days = Number_of_days_in_time_series-1;
    for (int i = 0; i <Number_of_days_in_time_series; i++){
        for (int j = 0; j < Number_of_variables_in_time_series; j++){
            DataPredict->Data_Time_Series[i][j] = Data[i][j];
        }
    }
    
    //Generation of initial individuals
    for (int i = 0; i < popsize; i++) {
        set_individuals(&population[i]);
        CoreModelVersusDataQuadraticError(&population[i], DataPredict, Data);
    }
    
    //Identification of the best individual in the origiginal population
    individual * fittestInActualPop = fittest_ind(population);
    printIndividual(fittestInActualPop);
    
    //INITIATION OF THE GENETIC ALGORITHM
    printf("\n\n--------GENETIC ALGORITHM-------\n\n");
    for (int iter=0; iter<1000; iter++) {
        //printf("\n- - - - -GENERATION %d- - - - -\n", iter);
        //Generation and Memory allocation for the children generation created
        individual * nextGen = NULL;
        if ((nextGen = (individual*) malloc(popsize*sizeof(individual))) == NULL) ExitError("could not allocate memory for the next generation", 1);
    
        //Genetic information modification
        for (int j = 0; j < popsize; j+=2) {
            //Tournament selection to decide 2 parents with lower fitness
            individual * Parent1 = TournamentSelection(population);
            individual * Parent2 = TournamentSelection(population);
            //Crossover step (in take just the full parameters without crossover in the middle of none
            OnePointCrossover(Parent1, Parent2, &nextGen[j], &nextGen[j+1]);
            // 2 mutations introduced in each chill (we add 2 as we want it very exploratory)
            Mutation(&nextGen[j]); Mutation(&nextGen[j]);
            Mutation(&nextGen[j+1]); Mutation(&nextGen[j+1]);
        }
        
        //Storage of the new generation and supresion of the old one
        individual * aux = population;
        population = nextGen;
        free(aux);
        
        //Calculation of the data for all individuals of the new generation and calculation of its fitness
        for (int i = 0; i < popsize; i++) CoreModelVersusDataQuadraticError(&population[i], DataPredict, Data);
        
        //Identification of the best individual in the new population
        fittestInActualPop = fittest_ind(population);
        //printIndividual(fittestInActualPop);
    }
    
    //Impesion of the best individual found in the last generation in order to compare it with the inital one
    printIndividual(fittestInActualPop);
    return 0;
    
    free(population);
    free(DataPredict);
};
