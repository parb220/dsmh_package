#include <fstream>
#include <cmath>

#include "prcsn.h"
#include "CSampleIDWeight.hpp"
#include "CEESParameter.hpp"

using namespace std; 

CEESParameter::CEESParameter() : 
	storage_dir(string()), 
	storage_marker(0), 
	run_id(string()), 
	number_energy_stage(0),
	number_striation(0), 
	pee(0.0),
	lambda_1(0.0), 
	lambda(vector<double>(0)),
	highest_stage(0), 
	lowest_stage(0), 
	THIN(0), 
	simulation_length(0), 
	burn_in_length(0)
{}

CEESParameter::~CEESParameter()
{}		

bool CEESParameter::SetTemperature_geometric()
// Given:
//	lambda[number_energy_stage] = 0.0; 
// 	lambda[number_energy_stage-1] = lambda_1
// 	lambda[0] = 1.0
//
// We first determine r = lambda[i]/lambda[i-1] for 1<=i <=number_energy_stage-1
// and then determine lambda for stages from 1 through number_energy_stage-2	
{
	// lambda
	if (number_energy_stage <= 0)
		return false; 
	lambda.resize(number_energy_stage+1); 
	lambda[0] = 1.0;
	lambda[number_energy_stage] = 0.0;  
	if (number_energy_stage > 1)
	{
		double r = log(lambda_1)/(number_energy_stage-1); 
		for (int i=1; i<number_energy_stage; i++)
			lambda[i] = exp((double)i*r); 	
	}
	return true; 
}

bool CEESParameter::SetTemperature_quadratic()
// Given:
// 	lambda[0] = 1.0
//	lambda[number_energy_stage] = 0.0; 
// Determine lambda as 
// 	lambda[i] = (H-i)^2/H^2
// Does not need information of lambda_1
{
	// lambda
	if (number_energy_stage <= 0)
		return false; 
	lambda.resize(number_energy_stage+1); 
	lambda[0] = 1.0;
	lambda[number_energy_stage] = 0.0; 
	for (int i=1; i<number_energy_stage; i++)
		lambda[i] = ((double)(number_energy_stage-i) * (double)(number_energy_stage-i)) /((double)number_energy_stage*(double)number_energy_stage ); 
	return true; 
}

bool CEESParameter::SetTemperature_polynomial(double r)
// Given:
// 	lambda[0] = 1.0
// 	lambda[number_energy_stage] = 0.0
// Determine lambda as
// 	lambda[i] = (H-i)^r/H^r
// where r needs to be specified by the user	
{
	if (number_energy_stage <= 0)
		return false; 	
	lambda.resize(number_energy_stage+1); 
	lambda[0] = 1.0; 
	lambda[number_energy_stage] = 0.0; 
	double denominator = pow((double)number_energy_stage,r); 
	for (int i=1; i<number_energy_stage; i++)
		lambda[i] = pow((double)(number_energy_stage-i),r)/denominator; 
	return true; 
}

double CEESParameter::LogRatio_Stage(const CSampleIDWeight &x, const CSampleIDWeight &y, int stage) const 
// Given two draws and the index of stage, calculate the LogRatio of posterior probabilities, where the likelihood 
// part of the posterior probability is raised to the power of lambda[stage]
// This log ratio is used to determine whether to accept an EE jump
{
	double log_prob_x_bounded = lambda[stage]*x.reserved + (x.weight-x.reserved); 
	double log_prob_y_bounded = lambda[stage]*y.reserved + (y.weight-y.reserved); 
        return log_prob_x_bounded - log_prob_y_bounded;
}


