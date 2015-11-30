#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <mpi.h>
#include <getopt.h>
#include <vector>
#include <iomanip>
#include "dw_rand.h"
#include "sbvar.hpp"
#include "CEquiEnergy_TimeSeries.hpp"
#include "CSampleIDWeight.hpp"
#include "CEESParameter.hpp"
#include "CStorageHead.hpp"
#include "mpi_constant.hpp"
#include "CMetropolis.hpp"
#include "storage_constant.hpp"
#include "TaskScheduling.hpp"
#include "maximization_option.hpp"

using namespace std; 

// bool WriteBlockScheme(const string &filename, const vector<TIndex> &blocks);

int main(int argc, char **argv)
{
	static struct option long_options[] =
        {
                {"data_file", required_argument, 0, 'd'},
                {"restriction_file", required_argument, 0, 's'},
		// Simulation options
		{"RunID", required_argument, 0, 'r'},
                {"Original", no_argument, 0, 'o'},
                {"Number of stages", required_argument, 0, 'n'},
		{"Number of striations", required_argument, 0, 'N'},
                {"Lambda_1", required_argument, 0, 'T'},
		{"Pee", required_argument, 0, 'P'}, 
                {"THIN", required_argument, 0, 'I'},
                {"ndraws", required_argument, 0, 'w'},
		{"burn-in length", required_argument, 0, 'W'}, 
                {"nInitial", required_argument, 0, 'S'},
                {"Highest Level", required_argument, 0, 'L'},
                {"Lowest Level", required_argument, 0, 'l'},
		// Diaoganistic options
		{"Number of groups for NSE analysis", required_argument, 0, 'G'}, 
		{0, 0, 0, 0}
	}; 

	int option_index = 0;
        size_t n_initial = 1;
	bool if_original = false;
	
	string data_file_name, restriction_file; 
	CEESParameter sim_option;
        sim_option.storage_marker = 10000;
        sim_option.run_id = to_string(time(NULL));
	sim_option.storage_dir = getenv("HOME")+string("/DW_TZ_GIT/projects_dw/work/sbvar/results/");
	// sim_option.storage_dir = getenv("HOME")+string("/work/mdd/results/");
        sim_option.lambda_1 = 0.1;
        sim_option.THIN = 50;
        sim_option.pee = 1.0/(10.0*sim_option.THIN);
        sim_option.simulation_length = 200000;
	sim_option.number_energy_stage = sim_option.number_striation = 1; 
	sim_option.highest_stage = sim_option.lowest_stage = -1; 

	// int block_scheme = 0; 

	int nGroup_NSE = 1; 

	while (1)
        {
                int c = getopt_long(argc, argv, "d:s:r:on:N:T:P:I:w:W:S:L:l:G:", long_options, &option_index);
                if (c == -1)
                        break;
		switch(c)
                {
                        case 'd':
                                data_file_name = string(optarg); break;
			case 's':
                                restriction_file = string(optarg); break;
			case 'r':
                                sim_option.run_id = string(optarg); break;
                        case 'o':
                                if_original = true; break; 
                        case 'n':
				sim_option.number_energy_stage = atoi(optarg); break; 
			case 'N':
				sim_option.number_striation = atoi(optarg); break; 
			case 'T':
                                sim_option.lambda_1 = atof(optarg); break;
			case 'P':
				sim_option.pee = atof(optarg); break; 
                        case 'I':
                                sim_option.THIN = atoi(optarg); break;
                        case 'w':
                                sim_option.simulation_length = atoi(optarg); break;
			case 'W': 
				sim_option.burn_in_length = atoi(optarg);  break; 
                        case 'S':
                                n_initial = atoi(optarg); break;
                        case 'L':
                                sim_option.highest_stage = atoi(optarg); break; 
                        case 'l':
                                sim_option.lowest_stage = atoi(optarg); break; 
			case 'G': 
				nGroup_NSE = atoi(optarg); break; 
			default: 
				break; 
		}
	}		
	if (data_file_name.empty() || restriction_file.empty()) 
	{
		cerr << "Usage: " << argv[0] << " -d data file -s restriction file.\n"; 
		abort(); 
	}	

	if (if_original)
        {
                sim_option.number_energy_stage = 1;
                sim_option.highest_stage = sim_option.lowest_stage = 0;
		sim_option.pee = 0.0; 
        }
	else 
	{
		if (sim_option.highest_stage < 0)
			sim_option.highest_stage = sim_option.number_energy_stage-1; 
		else 
        		sim_option.highest_stage = sim_option.highest_stage < sim_option.number_energy_stage-1 ? sim_option.highest_stage : sim_option.number_energy_stage-1; 
		if (sim_option.lowest_stage < 0)
			sim_option.lowest_stage = 0; 
		else 
        		sim_option.lowest_stage = sim_option.lowest_stage > 0 ? sim_option.lowest_stage : 0; 
	}
	// sim_option.SetTemperature_geometric(); // geometric 
	// else
	sim_option.SetTemperature_quadratic();  
	// sim_option.SetTemperature_polynomial(3.15); 

      	//////////////////////////////////////////////////////
      	// restrictions file
      	ifstream input; 
      	input.open(restriction_file.c_str(), ifstream::in);
      	if (!input.is_open()) 
	{
		cerr << "Unable to open restrictions file" << restriction_file << endl; 
		abort(); 
	}
      	vector<TDenseMatrix> U, V;
      	TDenseMatrix TrueA0, TrueAplus;
      	SetupRestrictionMatrices(U,V,TrueA0,TrueAplus,input);
      	input.close();
	
      	//////////////////////////////////////////////////////
      	// constant term, n_vars, n_predetermined, n_lags, n_exogenous, n_parameters
      	bool IsConstant=true;
      	int n_lags=NumberLags(TrueA0,TrueAplus,IsConstant); 
	int n_vars=TrueA0.cols; 
      	int n_parameters=0;
      	for (int i=0; i < n_vars; i++) 
		n_parameters+=U[i].cols+V[i].cols;

      	///////////////////////////////////////////////////////
      	//Generate/read data
      	int n_obs= 318-n_lags; // since 1988.01 
      	TDenseMatrix rawdata;

	input.open(data_file_name.c_str(), ifstream::in);
      	if (!input.is_open()) 
	{
		cout << "Unable to open data file" << data_file_name <<endl;
		abort(); 
	}
      	rawdata.Resize(n_obs+n_lags,1+n_vars);
      	input >> rawdata;
      	input.close();
	
      	TData_predetermined Data(n_lags, IsConstant, rawdata, TIndex(1,rawdata.cols-1), TIndex(), n_lags,rawdata.rows-1);
	
      	// Sims-Zha prior 
      	TDenseVector mu(6);
       	// mu(0)=1.0; mu(1)=1.0; mu(2)=1.0; mu(3)=1.2; mu(4)=3.0; mu(5)=3.0; 
       	mu(0) = 0.7; mu(1) = 0.5; mu(2) = 0.1; mu(3) = 1.2; mu(4) = 1.0; mu(5) = 1.0;
      	double periods_per_year=12.0;
      	SBVAR_symmetric_linear sbvar(&Data,mu,periods_per_year, 1.0, U, V);
      	sbvar.SetParameters(TrueA0,TrueAplus);

	// Estimate parameters
	sbvar.MaximizePosterior(1.0e-5,false);
     	TDenseVector EstimatedParameters(sbvar.NumberParameters());
      	sbvar.GetParameters(EstimatedParameters.vector);

	// blocking scheme
	/*vector <TIndex>blocks = sbvar.ConstructBlocks(block_scheme); 
	if (blocks.empty())
	{
		cerr << "block scheme can only be 0, 1 and 2.\n"; 
		abort(); 
	}*/

	/////////////////////////////////////////////////////////////////////
	// EquiEnergyModel
	MPI_Init(&argc, &argv);
        int my_rank, nNode;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nNode);
	
	dw_initialize_generator(time(NULL));

	CEquiEnergy_TimeSeries simulation_model;
        simulation_model.target_model = &sbvar;
        simulation_model.timer_when_started = -1;
        if (if_original)
                simulation_model.if_bounded = false;
        simulation_model.metropolis = new CMetropolis(&simulation_model); 
        simulation_model.parameter = &sim_option;
        simulation_model.current_sample = CSampleIDWeight(EstimatedParameters, 0, sbvar.LogPosterior(EstimatedParameters.vector), true);
        CSampleIDWeight mode = simulation_model.current_sample;
	simulation_model.storage = new CStorageHead (my_rank, sim_option.run_id, sim_option.storage_marker, sim_option.storage_dir, sim_option.number_energy_stage);

	// for communicating lower energy bound, and number of jumps from i-th striation to the j-th striation
	const int N_MESSAGE = (RESERVE_INDEX_START +1) + (sim_option.number_striation+1) + sim_option.number_striation*sim_option.number_striation; 

	if (my_rank == 0)
        {
		cout << "Lambda : " << endl; 
		for (int i=0; i<(int)(sim_option.lambda.size()); i++)
			cout << setprecision(20) << sim_option.lambda[i] << "\t"; 
		cout << endl; 
                if (!simulation_model.storage->makedir())
                {
                        cerr << "Error in making directory for " << sim_option.run_id << endl;
                        double *sMessage= new double [N_MESSAGE];
                        for (int i=1; i<nNode; i++)
                                MPI_Send(sMessage, N_MESSAGE, MPI_DOUBLE, i, END_TAG, MPI_COMM_WORLD);
                        delete [] sMessage;
                        exit(1);
                }
                
		// direct calculation of logMDD
		if (sim_option.highest_stage == sim_option.number_energy_stage-1)
		{
			for (int stage = sim_option.highest_stage; stage >= sim_option.lowest_stage; stage--)
			{
				sbvar.SetTemperature(sim_option.lambda[stage],1.0); 
				cout << "true logMDD at stage " << stage << setprecision(20) << ": " << sbvar.LogPosteriorIntegral() << endl; 
			}
			sbvar.SetTemperature(1.0,1.0); 
		}
		master_deploying(N_MESSAGE, nNode, n_initial, simulation_model, mode, nGroup_NSE); 
        }
        else
		 slave_computing(N_MESSAGE, simulation_model, mode); 

  	return 0;
}
