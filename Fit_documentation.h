#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "RKF78.h"
#include "RKF78.c"

//CONSTANTS USED
#define IC_GENES_NUMBER 3 //IC[] = { E(0), I{1}(0), A(0)}
#define PARAMETERS_GENES_NUMBER 11 //Pars = { beta, phi, epsI, epsY, sigma, gamma1, gamma2, kappa, p, alpha, delta}}
#define CoreModelDIM 8 //For Range-Kutta
#define Number_of_days_in_time_series 101
#define Number_of_variables_in_time_series 5
#define HMAX 1.0
#define HMIN 1.e-3
#define RKTOL 1.e-5

//OPERATIONS USED TO TRANSFORM PHENOTYPE IN GENOTYPE
#define crom2IC(c) (((double) c)/1000)
#define crom2HSPar(c) (((double) c)/1099511627776UL)
#define crom2Par(c) (((double) c)/1048576U)
#define crom2LSPar(c) (((double) c)/1024U)

//DEFINITION OF STRUCTURES USED
typedef struct {
    unsigned long IC[IC_GENES_NUMBER];
    unsigned long Pars[PARAMETERS_GENES_NUMBER];
    double fitness;
} individual;

typedef struct {
    double beta, phi, epsI, epsY;
    double sigma;
    double gamma1, gamma2;
    double kappa;
    double p;
    double alpha;
    double delta;
    double PopSize;
} ODE_Parameters;

typedef struct {
    double PopSize;
    unsigned N_days;
    double Data_Time_Series[Number_of_days_in_time_series][Number_of_variables_in_time_series];
} DataForFitting;

//Function where the derivates to be solved are expressed
void CoreModel(double t, double *x, unsigned dim, double *der, void * Params){
    ODE_Parameters * par = (ODE_Parameters *) Params;
    double sigmae = par->sigma*x[1], gamma1i1 = par -> gamma1*x[2], kappaA = par -> kappa*x[3], alphai2 = par -> alpha*x[5];
    der[0] = par->phi*x[2] + x[3] + (1-par->epsI)*(x[4]+x[5]) + (1-par -> epsY)*x[6];
    der[0] = -par->beta * (x[0] + der[0])/par->PopSize;
    der[1] = -der[0] - sigmae;
    der[2] = sigmae - gamma1i1;
    der[3] = (1-par->p)*gamma1i1 - kappaA - par->gamma2*x[3] ;
    der[4] = kappaA - par->gamma2*x[4];
    der[5] = par->p*gamma1i1 - par->gamma2*x[5] - alphai2;
    der[6] = alphai2 - (par->gamma2+par->delta)*x[6];
    der[7] = par->gamma2*(x[3] + x[4] + x[5] + x[6]);
}

//FITNES FUNCTION: Simple siumatory of the differences of all terms
double fitness_function(DataForFitting * Pred ,double Data[101][5], register unsigned ndays, double acc_error){
    acc_error += fabs(Pred->Data_Time_Series[ndays][0]- Data[ndays][0]) + fabs(Pred->Data_Time_Series[ndays][1]-Data[ndays][1]) + fabs(Pred->Data_Time_Series[ndays][2]-Data[ndays][2]) + fabs(Pred->Data_Time_Series[ndays][3]-Data[ndays][3]) + fabs(Pred->Data_Time_Series[ndays][4]-Data[ndays][4]);
    return acc_error;
}

//Function that uses Range_Kutta78 to calculate all the data for all day and all the EDO parametres for each individual (we use the function RKF78Sys which has been given by theory)
int GeneratePredictionFromIndividual(double *xt, void *ODE_pars, DataForFitting *Pred, double Data[101][5]) {
    register unsigned ndays; double fitness_error = 0;
    double t = 0.0, err, h = 1.e-3; int iter = 1;
    for(ndays=1; ndays <= Pred->N_days; ndays++) { int status;
        while(t+h < ndays) {
            status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
            if(status) return status;
            iter += 1;
        }
        
        h = ndays - t;
        status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);

        if(status) return status;

        Pred->Data_Time_Series[ndays][0] = xt[4];
        Pred->Data_Time_Series[ndays][1] = xt[5];
        Pred->Data_Time_Series[ndays][2] = xt[6];
        Pred->Data_Time_Series[ndays][3] = xt[7];
        Pred->Data_Time_Series[ndays][4] = Pred->PopSize - (xt[0]+xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]+xt[7]);
        fitness_error = fitness_function(Pred, Data, ndays, fitness_error);
    }
    return fitness_error;
}

//Fuction that stores the needed values for calculate all data in vectors. Then GeneratePredictionFromIndividual in introduced to do that calculation
void CoreModelVersusDataQuadraticError(individual *ind, void *TheData, double Data[101][5]) {
    double fitness_sum;
	DataForFitting *TDfF = (DataForFitting *) TheData;
	DataForFitting ThePrediction = { TDfF->PopSize, TDfF->N_days, {
		{ TDfF->Data_Time_Series[0][0], TDfF->Data_Time_Series[0][1],
		  TDfF->Data_Time_Series[0][2], TDfF->Data_Time_Series[0][3],
		  TDfF->Data_Time_Series[0][4] } } };

	double xt[CoreModelDIM] = { TDfF->PopSize, crom2IC(ind->IC[0]), crom2IC(ind->IC[1]),
		crom2IC(ind->IC[2]), TDfF->Data_Time_Series[0][0],
		TDfF->Data_Time_Series[0][1], TDfF->Data_Time_Series[0][2],
		TDfF->Data_Time_Series[0][3] };
		xt[0] -= (xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]);
    
	ODE_Parameters ODE_pars = { crom2HSPar(ind->Pars[0]), crom2Par(ind->Pars[1]),
			            crom2LSPar(ind->Pars[2]), crom2LSPar(ind->Pars[3]),
			            crom2Par(ind->Pars[4]), crom2Par(ind->Pars[5]),
				    crom2Par(ind->Pars[6]), crom2LSPar(ind->Pars[7]),
			            crom2Par(ind->Pars[8]), crom2HSPar(ind->Pars[9]),
				    crom2HSPar(ind->Pars[10]), TDfF->PopSize };
        
    fitness_sum = GeneratePredictionFromIndividual(xt, &ODE_pars, &ThePrediction, Data);
	ind->fitness = fitness_sum;
}

//ORIGINAL DATA OF THE EVOLUTION OF THE PANDEMIC
double Data[Number_of_days_in_time_series][Number_of_variables_in_time_series] = {
{1.000, 1.000, 0.000, 0.000, 0.000},
{1.841, 1.253, 0.056, 0.348, 0.000},
{2.285, 1.607, 0.123, 0.733, 0.002},
{2.571, 2.041, 0.203, 1.169, 0.004},
{2.812, 2.541, 0.300, 1.670, 0.008},
{3.059, 3.096, 0.414, 2.249, 0.013},
{3.337, 3.702, 0.546, 2.915, 0.020},
{3.655, 4.356, 0.698, 3.681, 0.029},
{4.018, 5.060, 0.870, 4.555, 0.040},
{4.428, 5.814, 1.062, 5.548, 0.054},
{4.889, 6.626, 1.276, 6.673, 0.071},
{5.402, 7.499, 1.513, 7.941, 0.091},
{5.972, 8.442, 1.773, 9.364, 0.115},
{6.603, 9.462, 2.058, 10.959, 0.143},
{7.299, 10.570, 2.369, 12.739, 0.175},
{8.067, 11.776, 2.710, 14.723, 0.212},
{8.914, 13.091, 3.083, 16.931, 0.253},
{9.847, 14.529, 3.490, 19.382, 0.301},
{10.875, 16.103, 3.934, 22.102, 0.355},
{12.008, 17.830, 4.420, 25.116, 0.415},
{13.257, 19.725, 4.952, 28.454, 0.483},
{14.633, 21.808, 5.534, 32.147, 0.559},
{16.150, 24.099, 6.172, 36.231, 0.644},
{17.823, 26.620, 6.871, 40.746, 0.738},
{19.667, 29.395, 7.638, 45.735, 0.844},
{21.701, 32.453, 8.480, 51.245, 0.960},
{23.943, 35.822, 9.406, 57.331, 1.090},
{26.416, 39.535, 10.422, 64.051, 1.233},
{29.143, 43.628, 11.540, 71.469, 1.392},
{32.151, 48.141, 12.770, 79.656, 1.569},
{35.469, 53.116, 14.124, 88.693, 1.763},
{39.128, 58.603, 15.615, 98.666, 1.979},
{43.164, 64.654, 17.256, 109.670, 2.217},
{47.616, 71.328, 19.064, 121.812, 2.480},
{52.527, 78.688, 21.056, 135.210, 2.770},
{57.944, 86.805, 23.252, 149.991, 3.091},
{63.918, 95.758, 25.671, 166.298, 3.446},
{70.508, 105.632, 28.338, 184.289, 3.837},
{77.777, 116.523, 31.278, 204.136, 4.269},
{85.795, 128.535, 34.520, 226.031, 4.745},
{94.638, 141.782, 38.094, 250.183, 5.271},
{104.391, 156.393, 42.035, 276.826, 5.851},
{115.149, 172.506, 46.380, 306.216, 6.492},
{127.013, 190.277, 51.172, 338.635, 7.198},
{140.099, 209.874, 56.456, 374.394, 7.978},
{154.529, 231.485, 62.282, 413.837, 8.838},
{170.444, 255.316, 68.708, 457.342, 9.787},
{187.994, 281.594, 75.793, 505.327, 10.833},
{207.347, 310.569, 83.606, 558.252, 11.988},
{228.687, 342.516, 92.222, 616.623, 13.261},
{252.217, 377.737, 101.722, 681.000, 14.66},
{278.160, 416.566, 112.197, 751.997, 16.215},
{306.763, 459.370, 123.746, 830.294, 17.924},
{338.296, 506.551, 136.480, 916.637, 19.809},
{373.057, 558.552, 150.519, 1011.850, 21.887},
{411.372, 615.862, 165.995, 1116.839, 24.180},
{453.603, 679.015, 183.056, 1232.602, 26.708},
{500.144, 748.598, 201.861, 1360.238, 29.496},
{551.431, 825.258, 222.588, 1500.956, 32.570},
{607.940, 909.702, 245.431, 1656.087, 35.960},
{670.196, 1002.705, 270.602, 1827.096, 39.697},
{738.773, 1105.119, 298.338, 2015.593, 43.818},
{814.301, 1217.874, 328.894, 2223.347, 48.361},
{897.473, 1341.989, 362.553, 2452.305, 53.369},
{989.042, 1478.579, 399.625, 2704.605, 58.889},
{1089.838, 1628.858, 440.449, 2982.593, 64.974},
{1200.765, 1794.155, 485.396, 3288.848, 71.680},
{1322.811, 1975.914, 534.874, 3626.195, 79.069},
{1457.053, 2175.710, 589.326, 3997.735, 87.212},
{1604.666, 2395.253, 649.239, 4406.864, 96.183},
{1766.929, 2636.396, 715.142, 4857.302, 106.065},
{1945.231, 2901.150, 787.612, 5353.119, 116.949},
{2141.079, 3191.686, 867.280, 5898.765, 128.936},
{2356.106, 3510.345, 954.827, 6499.101, 142.133},
{2592.078, 3859.646, 1050.997, 7159.429, 156.661},
{2850.898, 4242.291, 1156.593, 7885.531, 172.651},
{3134.614, 4661.169, 1272.484, 8683.701, 190.246},
{3445.426, 5119.362, 1399.610, 9560.781, 209.600},
{3785.685, 5620.137, 1538.981, 10524.199, 230.885},
{4157.899, 6166.949, 1691.682, 11582.004, 254.286},
{4564.733, 6763.432, 1858.876, 12742.907, 280.004},
{5009.003, 7413.384, 2041.802, 14016.312, 308.258},
{5493.676, 8120.751, 2241.781, 15412.350, 339.286},
{6021.857, 8889.602, 2460.210, 16941.914, 373.345},
{6596.777, 9724.093, 2698.565, 18616.677, 410.714},
{7221.773, 10628.429, 2958.390, 20449.122, 451.691},
{7900.265, 11606.813, 3241.302, 22452.549, 496.601},
{8635.723, 12663.380, 3548.971, 24641.082, 545.789},
{9431.630, 13802.127, 3883.121, 27029.662, 599.628},
{10291.435, 15026.829, 4245.506, 29634.031, 658.513},
{11218.496, 16340.939, 4637.900, 32470.693, 722.867},
{12216.019, 17747.481, 5062.074, 35556.871, 793.138},
{13286.985, 19248.930, 5519.771, 38910.428, 869.799},
{14434.065, 20847.086, 6012.679, 42549.783, 953.349},
{15659.535, 22542.936, 6542.399, 46493.794, 1044.309},
{16965.181, 24336.518, 7110.407, 50761.622, 1143.224},
{18352.193, 26226.782, 7718.015, 55372.566, 1250.660},
{19821.068, 28211.456, 8366.327, 60345.878, 1367.197},
{21371.503, 30286.926, 9056.199, 65700.547, 1493.434},
{23002.296, 32448.131, 9788.183, 71455.071, 1629.977},
{24711.260, 34688.480, 10562.489, 77627.195, 1777.4}};


