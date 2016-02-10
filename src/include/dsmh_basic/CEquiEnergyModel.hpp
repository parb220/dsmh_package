#ifndef CLASS_EQUI_ENERGY_MODEL_HEADER
#define CLASS_EQUI_ENERGY_MODEL_HEADER

#include <vector>
#include <ctime>
#include <cstdio>
#include "CSampleIDWeight.hpp"
#include "CEESParameter.hpp"
#include "CStorageHead.hpp"

using namespace std;

const int METROPOLIS_JUMP = 1; 
const int EQUI_ENERGY_JUMP = -1; 
const int NO_JUMP = 0; 

class CSampleIDWeight;
class CStorageHead;
class CEESParameter;
class CMetropolis; 
class CEquiEnergyModel; 
class TDenseVector;
class TDenseMatrix;


class CEquiEnergyModel
{
private:
	CEquiEnergyModel(const CEquiEnergyModel &);
	std::vector<TDenseVector> gmm_mean; 
	std::vector<TDenseMatrix> gmm_covariance_sqrt; 
	std::vector<double> gmm_covariance_sqrt_log_determinant; 
	std::vector<TDenseMatrix> gmm_covariance_sqrt_inverse; 

public:
///////////////////////////////////////////////////////////////
// Parameters for equi-energy sampling
	bool if_bounded;	// whether likelihood should be raised to power of lambda 
 	int energy_stage;	// index of current stage	
	double lambda; 		// lambda of current stage

//////////////////////////////////////////////////////////////
// Current sample 
	CSampleIDWeight current_sample; 
	int timer_when_started; 
	
public:
	CMetropolis *metropolis; 	// pointer to CMetropolis
	CEESParameter *parameter;	// pointer to CEESParameter 
	CStorageHead *storage; 		// pointer to CStorageHead
	
	virtual double log_posterior_function(CSampleIDWeight &x)=0;
        virtual double log_likelihood_function(const CSampleIDWeight &x)=0;
	virtual double log_prior_function(const CSampleIDWeight &x) = 0; 
	virtual double log_posterior_function(const double *x, int n)=0; 
        virtual double log_likelihood_function(const double *x, int n)=0;
	virtual double log_prior_function(const double *x, int n) = 0;  

	// Draw samples
	int EE_Draw(); 	// equi-energy draw
	bool MakeEquiEnergyJump(CSampleIDWeight &y_end, const CSampleIDWeight &y_initial); 
	bool JumpAcrossStriation(const CSampleIDWeight &y_end, const CSampleIDWeight &y_initial, TDenseMatrix &) const; 
	virtual void SaveSampleToStorage(const CSampleIDWeight &sample); 
	virtual void Take_New_Sample_As_Current_Sample(const CSampleIDWeight &x); 

public:
	std::vector<int> BurnIn(int burn_in_length) ;	// Return a two-element vector, containing the number of EE jumps and the number of MH jumps
	std::vector<int> Simulation_Within(TDenseMatrix &jump_table, bool if_storage, const string &sample_file_name=string()); 	// Simulation within the same energy stage (no jumping across stages). Return a two-element vector, containing the number of EE jumps and the number of MH jumps
	std::vector<int> Simulation_Cross(TDenseMatrix &jump_table, bool if_storage, const string &sample_file_name=string()); 	// Simulation across stages. Return a two-element vector, containing the number of EE jumps and the number of MH jumps	
	std::vector<int> Simulation_Prior(bool if_storage, const string &sample_file_name=string()); 	// Simulation from prior. Return a two-element vector, containing the number of EE jumps and the number of MH jumps 
	
	virtual bool DrawParametersFromPrior(double *x) const = 0;
	
	std::vector<CSampleIDWeight> Initialize_WeightedSampling(const std::vector<CSampleIDWeight> &, int K, int stage_index) ; // Make K draws from the draws of the previous stage (as the input) according to the re-calculated weight. Return the K draws
	// Reweight samples
	vector<double> Reweight(const vector<CSampleIDWeight> &samples, int current_stage, int previous_stage); // Recalculate of the samples according to the lambda's of the current and prevsious stage indices

///////////////////////////////////////////////////////////////////////////////////////////
// Construction & Destruction functions here
public:
	CEquiEnergyModel(); 
	CEquiEnergyModel(bool _if_bounded, int eL, double _t, const CSampleIDWeight & _x=CSampleIDWeight(), time_t _time=time(NULL), CMetropolis *_metropolis =NULL, CEESParameter *_parameter=NULL, CStorageHead *_storage = NULL); 
	virtual ~CEquiEnergyModel() {} 
};

#endif
