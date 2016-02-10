#ifndef CLASS_EES_PARAMETER_HEADER
#define CLASS_EES_PARAMETER_HEADER

#include <string>
#include <vector>

using namespace std; 

class CSampleIDWeight; 

class CEESParameter 
{
public:
	string storage_dir;	// fold to store 
	int storage_marker;	// number of samples stored in memory before replenish from the hard drive
	string run_id;		// id of the run 
	int number_energy_stage;	// number of energy stages
	int number_striation; 	// number of rings
	double pee; 	// probability of equi-energy-jump
	double lambda_1;	// lambda for the highest temperature stage
	vector<double> lambda; 	// lambda for all temperature stages
	
	int highest_stage; 
	int lowest_stage; 
	int THIN;  	// thinning factor for MH samples
	int simulation_length; // length of simulation
	int burn_in_length; // length of burn-in
public:
	CEESParameter(); 
	~CEESParameter(); 
	
	bool SetTemperature_geometric(); 
	bool SetTemperature_quadratic();  
	bool SetTemperature_polynomial(double r);

	double LogRatio_Stage(const CSampleIDWeight &x, const CSampleIDWeight &y, int stage) const; 
};
#endif
