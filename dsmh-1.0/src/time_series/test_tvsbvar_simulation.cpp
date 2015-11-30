#include <cstdlib>
#include <iomanip>
#include <getopt.h>
#include <fstream>
#include <ctime>
#include <mpi.h>
#include "dw_rand.h"
#include "dw_dense_matrix.hpp"
#include "dw_data.hpp"
#include "regime_processes.hpp"
#include "tvsbvar.hpp"
#include "CEESParameter.hpp"
#include "CSampleIDWeight.hpp"
#include "CStorageHead.hpp"
// #include "maximization_option.hpp"
#include "mpi_constant.hpp"	
#include "CEquiEnergy_TimeSeries.hpp"
#include "CMetropolis.hpp"
#include "storage_constant.hpp"
#include "TaskScheduling.hpp"

using namespace std; 

string cluster_to_string(int);

// bool WriteBlockScheme(const string &filename, const vector<TIndex> &blocks); 

int main(int argc, char **argv)
{
	// Initialize MPI
	MPI_Init(&argc, &argv);
        int my_rank, nNode;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nNode);

	// random number generator
	dw_initialize_generator(time(NULL)+my_rank*100000); 

	// command line options
	static struct option long_options[] =
        {
                {"data_file", required_argument, 0, 'd'},
		{"markov_regime_process_file", required_argument, 0, 'm'},
                {"deterministic_regime_process_file", required_argument, 0, 'c'},
		{"restriction_file", required_argument, 0, 's'},
		{"hyper_parameter_file", required_argument, 0, 'h'},
		{"initial_parameter_file", required_argument, 0, 'p'},
		// Simulation options
		{"RunID", required_argument, 0, 'r'},
		{"Original", no_argument, 0, 'o'},
		{"Number of stages", required_argument, 0, 'n'},
		{"Number of striations", required_argument, 0, 'N'},
		{"Lambda_1", required_argument, 0, 'T'},
		{"Pee", required_argument, 0, 'P'},
                {"THIN", required_argument, 0, 'I'},
                {"ndraws", required_argument, 0, 'w'},
		{"burn in Length", required_argument, 0, 'W'}, 
                {"nInitial", required_argument, 0, 'S'}, 
		{"Highest Level", required_argument, 0, 'L'}, 
		{"Lowest Level", required_argument, 0, 'l'},
		{"Number of groups for NSE analysis", required_argument, 0, 'G'},
		// HillClimb options
		/*{"block_scheme", required_argument, 0, 'b'}, 
		{"perturbation_scheme", required_argument, 0, 'e'}, 
		{"perturbation_iteration", required_argument, 0, 'E'}, 
		{"block_iteration", required_argument, 0, 'K'}, 
		{"maximization_iteration", required_argument, 0, 'M'},*/
                {0, 0, 0, 0}
        };

	int option_index = 0;
	int n_initial = 1; 
	bool if_original = false; 
	
	string data_file_name, markov_regime_process_file, deterministic_regime_process_file, restriction_file, hyper_parameter_file, initial_parameter_file; 
	
	CEESParameter sim_option;
	sim_option.storage_marker = 10000; 
	sim_option.run_id = cluster_to_string(time(NULL)); 
	sim_option.storage_dir = getenv("HOME")+string("/DW_TZ_GIT/projects_dw/work/tvsbvar/results/"); 
	sim_option.lambda_1 = 0.1;
	sim_option.pee = 1.0/(10.0*sim_option.THIN);
	sim_option.THIN = 50; 
	sim_option.simulation_length = 200000; // SIMULATION_LENGTH; 
	sim_option.number_energy_stage = sim_option.number_striation = 1; 
	sim_option.highest_stage = sim_option.lowest_stage = -1; 

	/*MaximizationOptions max_option; 
	max_option.BlockScheme = 0; 
	max_option.PerturbationScale = 1.0; 
	max_option.MaxPerturbationIterations = 10; 
	max_option.MaxBlockIterations = 30; 
	max_option.MaxOptimizationIterations = 10; 
	max_option.ConstantOptimization = false;  */

	int nGroup_NSE = 1;
	
	while (1)
        {
                int c = getopt_long(argc, argv, "d:m:c:s:h:p:r:on:N:t:T:P:I:w:W:S:L:l:G:b:e:E:K:M:O", long_options, &option_index);
                if (c == -1)
                        break;
                switch(c)
                {
                        case 'd':
                                data_file_name = string(optarg); break;
			case 'm': 
				markov_regime_process_file = string(optarg); break; 
			case 'c':
				deterministic_regime_process_file = string(optarg); break;
			case 's': 
				restriction_file = string(optarg); break; 
			case 'h':
				hyper_parameter_file = string(optarg); break; 
			case 'p': 
				initial_parameter_file = string(optarg); break; 
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
				sim_option.burn_in_length = atoi(optarg); break; 
			case 'S': 
				n_initial = atoi(optarg); break; 
			case 'L': 
				sim_option.highest_stage = atoi(optarg);  break; 
			case 'l': 
				sim_option.lowest_stage = atoi(optarg); break; 
			case 'G':
                                nGroup_NSE = atoi(optarg); break;
			/*case 'b':
				max_option.BlockScheme = atoi(optarg); break; 
			case 'e': 
				max_option.PerturbationScale = atof(optarg); break; 
			case 'E': 
				max_option.MaxPerturbationIterations = atoi(optarg); break; 
			case 'K':
				max_option.MaxBlockIterations = atoi(optarg); break; 
			case 'M':
				max_option.MaxOptimizationIterations = atoi(optarg); break; 
			case 'O': 
				max_option.ConstantOptimization = true; break;  */
			default:
                                break;
                }
        }
        if (data_file_name.empty() || (markov_regime_process_file.empty() && deterministic_regime_process_file.empty() )|| restriction_file.empty() || hyper_parameter_file.empty())
	{
		cerr << "Usage: " << argv[0] << " -d data file -m markov regime process file -c deterministic regime process file -s restriction file -h hyper parameter file\n"; 
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
        sim_option.SetTemperature_quadratic();

	/////////////////////////////////// Data //////////////////////////////////////////
	TDenseMatrix data_from_file(318,4);  //570,4); // data_from_file(90, 37); // 
        ifstream input_file;
        input_file.open(data_file_name.c_str(), fstream::in);
        if (!input_file)
        {
                cerr << "Error in opening " << data_file_name << endl;
                abort();
        }
        input_file >> data_from_file;
        input_file.close();
	// append a constant vector Ones(data_from_file.rows) to the end of raw_data
	TDenseMatrix raw_data(data_from_file.rows, data_from_file.cols+1);
        raw_data.Insert(0,0,data_from_file);
        raw_data.InsertColumnMatrix(0,data_from_file.cols,TDenseVector(data_from_file.rows,1.0));
	// Generate TData_predetermined from raw_data
	int n_lags = 13; // 5;
        TData_predetermined input_data(n_lags,false,raw_data, TIndex(1,data_from_file.cols-1), TIndex(data_from_file.cols), n_lags, data_from_file.rows-1);// TIndex(8)(2)(35)(12)(20)(1), TIndex(data_from_file.cols), n_lags+16, data_from_file.rows-1); //

	//////////////////////////// Regime process: independent markov ///////////////////////
	// markov
	TRegimeProcessArray regime_process_array;
        if (!markov_regime_process_file.empty() )
        {
                try { regime_process_array = SetupMarkovArrayRegimeProcess(markov_regime_process_file.c_str()); }
                catch (dw_exception &e)
                {
                        cerr << "Error in reading or parsing " << markov_regime_process_file << endl << e.what() << endl;
                        abort();
                }
        }
	// deterministic
	if (!deterministic_regime_process_file.empty())
        {
                try
                {
                        int rows=74;
                        input_file.open(deterministic_regime_process_file.c_str());
                        if (!input_file.is_open())
                                throw dw_exception("unable to open regimes file");
                        TDenseMatrix regime_file(rows,3);
                        input_file >> regime_file;
                        input_file.close();
                        vector<int> regimes(rows);
                        for (int i=0; i < 2; i++)
                        {
                                int n_regimes = 1;
                                for (int t=0; t<rows; t++)
                                {
                                        regimes[t]=regime_file(t,i)-1;
                                        if (regimes[t] < 0)
                                                throw dw_exception("regimes cannot be negative");
                                        if (regimes[t] >= n_regimes)
                                                n_regimes=regimes[t]+1;
                                }
                                regime_process_array.push_back(new TRegimeProcess_deterministic(n_regimes,regimes));
                        }
                }
                catch (dw_exception &e)
                {
                        cerr << "Error in reading or parsing " << deterministic_regime_process_file << endl << e.what() << endl;
                        abort();
                }
        }
	TRegimeProcess_independent independent_regime_process(regime_process_array);
        TDenseVector parameters(independent_regime_process.NumberParameters(),0.0);
        independent_regime_process.SetParameters(parameters.vector);

	/////////////////////////// Hyper parameter ////////////////////////////////////////
	input_file.open(hyper_parameter_file.c_str(), fstream::in);
        if (!input_file)
        {
                cerr << "Error in opening " << hyper_parameter_file << endl;
                abort();
        }
        THyperParameter_Vector hyper_parameter(input_file);
        input_file.close();
	
	//////////////////////////// tvsbvar: /////////////////////////////////////////////////
	input_file.open(restriction_file.c_str(), fstream::in);
        if (!input_file)
        {
                cerr << "Error in opening" << restriction_file << endl;
                abort();
        }

        Sims_Zha::log2pi = 1.837877066409345;
        TVSBVAR tvsbvar_model(input_file, (TData_predetermined *)(&input_data), &independent_regime_process, &hyper_parameter, 4, 1.0);
        input_file.close();

        TDenseVector parameters_initial(tvsbvar_model.NumberParameters(),0.0), parameters_optimal;
        double log_posterior_optimal;
        /*if (max_option.ConstantOptimization)
        {
                tvsbvar_model.Optimize_ConstantRegime_NPSOL(log_posterior_optimal, parameters_optimal, parameters_initial, max_option.PerturbationScale, max_option.MaxPerturbationIterations, 1);
                parameters_initial.CopyContent(parameters_optimal);
        }
        else
        {*/
                input_file.open(initial_parameter_file.c_str());
                if (input_file)
                {
                        input_file >> parameters_initial;
                        input_file.close();
                        tvsbvar_model.SetParameters(parameters_initial.vector);
                }
                else
                {
                        tvsbvar_model.DefaultParameters();
                        tvsbvar_model.GetParameters(parameters_initial.vector);
                }
         /*} */

	/////////////////////////////// max_option /////////////////////////////////////////////
	/*max_option.blocks = tvsbvar_model.ConstructBlocks(max_option.BlockScheme);
        if (max_option.blocks.empty())
        {
                cerr << "block scheme can only be 0, 1, 2, 3 or 4.\n";
                abort();
        }*/

	////////////////////////////  EquiEnergyModel ///////////////////////////////////////
	CEquiEnergy_TimeSeries simulation_model;
        simulation_model.target_model = &tvsbvar_model;
        simulation_model.timer_when_started = time(NULL);
        if (if_original)
                simulation_model.if_bounded = false;
        simulation_model.metropolis = new CMetropolis(&simulation_model); 
        simulation_model.parameter = &sim_option;
        simulation_model.current_sample = CSampleIDWeight(parameters_initial, 0, tvsbvar_model.LogPosterior(parameters_initial), true);
        CSampleIDWeight mode = simulation_model.current_sample;
	simulation_model.storage = new CStorageHead(my_rank, sim_option.run_id, sim_option.storage_marker, sim_option.storage_dir, sim_option.number_energy_stage); 

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
		master_deploying(N_MESSAGE, nNode, n_initial, simulation_model, mode, nGroup_NSE);
	}
	else
		// slave_computing(period, max_period, n_initial, simulation_model, mode, max_option.MaxOptimizationIterations, max_option.MaxPerturbationIterations, max_option.PerturbationScale);
		slave_computing(N_MESSAGE, simulation_model, mode); 
	return 0;
}
