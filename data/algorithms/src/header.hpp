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

#include <iostream>
#include <sstream>
#include <time.h>
#include <string.h>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <algorithm>
#include <iterator>   
#include <random>
#include <set>
#include <limits>
#include <map>
#include <iomanip>

using namespace std;

extern int num_pat, num_prv, num_prv_type, num_ser, num_age, num_pln, num_prf, num_zip, epsilon, sh_err;

struct profile {

	int age;
	int tob;
	int price_out;
	int pop = 0;

	double zipx;
	double zipy;

	double beta_prc;
	double beta_brd;

	vector<double> beta_qual;
	vector<double> beta_dist;
	vector<double> beta_wait;

	vector<double> pref;
	
	vector<int> glob_ch;
	vector<double> util_prv;
	vector<vector<double>*> util_prv_ser;
	vector<vector<double>*> pat_cost;
	vector<vector<double>*> pat_cost_self; 

};


struct provider {

	int type;
	int num_ser = 0;

	double zipx;
	double zipy;

	double qual;
	double ann_cost = 0;
	
	vector<double> wait;

	vector<double> neg_coeff;
	vector<double> self_coeff;

	vector<int> services;

};


struct plan {

	double single_profit = 0;

	vector<int> providers;
	vector<double> single_prem;
	vector<double> prem;
	vector<double> prem_comp;
	vector<double> cost;
	vector<double> self_cost;
	vector<double> brd_util;
	vector<double> util_comp;
	vector<double> util;

	plan() {
		single_prem = vector<double>(num_prf, 0);
		prem = vector<double>(num_prf, 0);
		prem_comp = vector<double>(num_prf, 0);
		cost = vector<double>(num_prf, 0);
		self_cost = vector<double>(num_prf, 0);
		util = vector<double>(num_prf, 0);
		util_comp = vector<double>(num_prf, 0);
		brd_util = vector<double>(num_prf, 0);
	}

};

extern vector<provider*> providers;
extern vector<plan*> plans;
extern vector<profile*> profiles;

extern string path_r, path_f;

void write_profiles_raw();
void write_profiles_fin();
void write_providers_raw();
void write_providers_fin();
void write_plans_fin();

bool load_data(string inst, vector<vector<int>*>& dat);
bool load_data(string inst, vector<vector<double>*>& dat);
bool load_data(string inst, vector<double>& dat);
bool load_data(string inst, vector<int>& dat);

int det_prof(int age, int distr, int tob);
double eval_plan(vector<int>& prvs);
void ord_providers(vector<int>& ord_prv);
void ord_providers_cost(vector<int>& ord_prv, bool low_to_high);
double round_p(double num, int prec);
string print_com(double value);
double calc_dist(double x1, double y1, double x2, double y2);


