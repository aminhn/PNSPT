/*
---------------------------------------------------------------------------------------------------------------------------------- 
Provider Network Selection and Patient Targeting in Health Insurance Markets
Copyright (C) 2022 Amin Hosseininasab
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
----------------------------------------------------------------------------------------------------------------------------------
*/

#pragma once
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

extern int num_pat, num_prv, num_pln, num_plnc, num_ser, tot_num_pat, epsilon, max_single_pln;
extern double max_single_prof, deduct;

extern string instance, in_path;


struct patient {
	
	int pop;
	int price_out;

	double beta_prc;
	double beta_brd;

	vector<int> glob_ch;
	vector<double> pref;
	vector<double> util_prv;
	vector<vector<double>*> util_prv_ser;
	vector<vector<double>*> cost_prv;
	vector<vector<double>*> cost_prv_self;

};


struct provider {

	long int ann_cost;
	int num_ser;
	vector<bool> services;

};


struct plan {

	vector<int> providers;
	vector<double> cost;
	vector<double> self_cost;
	vector<double> prem;
	vector<double> single_prem;
	vector<double> brd_util;
	vector<double> util;

};


extern vector<patient*> patients;
extern vector<provider*> providers;
extern vector<plan*> plans;
extern vector<plan*> comp_plans;

bool load_inst();

void calc_cost_util(plan* pln);



