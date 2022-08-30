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

#include "header.hpp"
#include "model.hpp"
#include <map>

using namespace std;

int num_pat, num_prv, num_prv_type, num_ser, num_age, num_pln, num_prf, num_zip, epsilon = 1800, ver = 0;
double tot_num_ser = 0, max_wait = -1, util_low = 10, util_high = 100000, max_dist = -1, c_anncost = 1, c_neg = 1, c_sneg = 1, c_cost = 1, c_dist = 1, c_pop = 1, c_age = 1, c_tob = 1, c_qual = 1, c_wait = 1, c_eps = 1, c_beta_pln = 1, c_beta_qual = 1, c_beta_dist = 1, c_beta_wait = 1;

bool same_ann = 0, same_cost = 0;

string path_r, path_f;

random_device rrand;
default_random_engine generator(rrand());

discrete_distribution<int> loc_dis;

vector<int> ser_import;				//1:important, 0:not
vector<double> pat_age_prob;			//0-18, 19-24, 25-44, 45-54, 55-64
vector<double> prv_type_prob;			//Amb Sur, Lab, Hospice, Clinic, Hospital, Office
vector<double> num_zip_prb;			//intercept, slope(per 100000)
vector<double> prv_qual_prob;			//low, med, high
vector<vector<double>*> prv_ser_prob;
vector<vector<double>*> pat_pref;		//rows: age_group, cols: ser_type - Hospice, Mat, Sur, PC EV M, ER, Lab
vector<vector<double>*> cost_ser;		//rows and cols as above
vector<vector<double>*> beta_pln_l;		//rows and cols as above
vector<vector<double>*> beta_prv_nl;		//rows and cols as above
vector<vector<double>*> beta_prv_sl;		//rows and cols as above
vector<vector<double>*> zip_loc;		//zipx,zipy,zipx std,zipy std
vector<vector<double>*> zip_denst;		//prb, std
vector<vector<double>*> tob_prb;

vector<double> zipx_loc;
vector<double> zipy_loc;

vector<vector<int>> prf_dem_pln;
vector<int> tot_dem_prv;

vector<provider*> providers;
vector<plan*> plans;
vector<profile*> profiles;

map<vector<int>,int> pat_to_prf;

void create_districts();
void create_profiles();
void create_patients();
void create_providers();
void calc_prf_util();
void create_plans();
void excl_heur();
void price_plan0();
void price_plans();
void incl_heur(plan* pln);


int main(int argc, char* argv[]) {

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] != '-')
			continue;
		else if (strcmp(argv[i], "-npat") == 0)
			num_pat = stoi(argv[i + 1]);
		else if (strcmp(argv[i], "-ver") == 0) {
			ver = stoi(argv[i + 1]);
			generator = default_random_engine(ver);
		}
		else if (strcmp(argv[i], "-c_ann") == 0)
			c_anncost = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_neg") == 0)
			c_neg = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_sneg") == 0)
			c_sneg = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_cst") == 0)
			c_cost = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_dist") == 0)
			c_dist = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_pop") == 0)
			c_pop = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_age") == 0)
			c_age = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_tob") == 0)
			c_tob = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_qual") == 0)
			c_qual = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_wait") == 0)
			c_wait = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_eps") == 0)
			c_eps = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_beta_pln") == 0)
			c_beta_pln = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_beta_qual") == 0)
			c_beta_qual = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_beta_dist") == 0)
			c_beta_dist = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-c_beta_wait") == 0)
			c_beta_wait = stod(argv[i + 1]);
		else if (strcmp(argv[i], "-s_ann") == 0)
			same_ann = 1;
		else if (strcmp(argv[i], "-s_cost") == 0)
			same_cost = 1;
		else if (strcmp(argv[i], "-s_all") == 0) {
			same_cost = 1;
			same_ann = 1;
		}


	}
	

	string mods = "";
	if (c_anncost != 1) {
		int temp = c_anncost * 10;
		mods = "ann" + to_string(temp) + "_";
	}
	if (c_neg != 1) {
		int temp = c_neg * 10;
		mods = "neg" + to_string(temp) + "_";
	}
	if (c_sneg != 1) {
		int temp = c_sneg * 10;
		mods = "slfneg" + to_string(temp) + "_";
	}
	if (c_cost != 1) {
		int temp = c_cost * 10;
		mods = "cost" + to_string(temp) + "_";
	}
	if (c_dist != 1) {
		int temp = c_dist * 10;
		mods = "dist" + to_string(temp) + "_";
	}
	if (c_wait != 1) {
		int temp = c_wait * 10;
		mods = "wait" + to_string(temp) + "_";
	}
	if (c_eps != 1) {
		int temp = c_eps * 10;
		mods = "eps" + to_string(temp) + "_";
		epsilon *= c_eps;
	}
	if (c_pop != 1) {
		int temp = c_pop * 10;
		mods = "pop" + to_string(temp) + "_";
	}
	if (c_age != 1) {
		int temp = c_age * 10;
		mods = "age" + to_string(temp) + "_";
	}
	if (c_tob != 1) {
		int temp = c_tob * 10;
		mods = "tob" + to_string(temp) + "_";
	}
	if (c_qual != 1) {
		int temp = c_qual * 10;
		mods = "qual" + to_string(temp) + "_";
	}
	if (c_beta_pln != 1) {
		int temp = c_beta_pln * 10;
		mods = "beta_pln" + to_string(temp) + "_";
	}
	if (c_beta_qual != 1) {
		int temp = c_beta_qual * 10;
		mods = "beta_qual" + to_string(temp) + "_";
	}
	if (c_beta_wait != 1) {
		int temp = c_beta_wait * 10;
		mods = "beta_wait" + to_string(temp) + "_";
	}
	if (c_beta_dist != 1) {
		int temp = c_beta_dist * 10;
		mods = "beta_dist" + to_string(temp) + "_";
	}
	if (same_ann == 1) 
		mods = "same_ann_";
	if (same_cost == 1) 
		mods = "same_cost_";
	if (same_cost == 1 && same_ann == 1) 
		mods = "same_all_";

	if (num_prv == 0)
		num_prv = round(238751.0 / 332403650 * num_pat);

	path_r = "./instances/" + mods + "I" + to_string(num_pat) + "_ver" + to_string(ver) + "/rawfiles";
	path_f = path_r.substr(0, path_r.length() - 8);
	string sys_exc1= "rm -rf " + path_f;
	string sys_exc2= "mkdir -p " + path_r;
	bool exc = system(sys_exc1.c_str());
	exc = system(sys_exc2.c_str());

	path_r = path_r + "/";

	string path_d = "./data/";
	
	load_data(path_d + "age.csv", pat_age_prob);
	load_data(path_d + "pat_pref.csv", pat_pref);
	load_data(path_d + "pat_costs.csv", cost_ser);
	load_data(path_d + "prv_type.csv", prv_type_prob);
	load_data(path_d + "prv_quality.csv", prv_qual_prob);
	load_data(path_d + "prv_serv.csv", prv_ser_prob);
	load_data(path_d + "impr_serv.csv", ser_import);
	load_data(path_d + "beta_prv_n.csv", beta_prv_nl);
	load_data(path_d + "beta_prv_s.csv", beta_prv_sl);
	load_data(path_d + "beta_pln.csv", beta_pln_l);
	load_data(path_d + "zip_loc.csv", zip_loc);
	load_data(path_d + "num_zip.csv", num_zip_prb);
	load_data(path_d + "zip_density.csv", zip_denst);
	load_data(path_d + "tobacco.csv", tob_prb);

	num_age = pat_age_prob.size();
	num_ser = cost_ser[0]->size();
	num_zip = round(num_zip_prb[0] + num_zip_prb[1] * num_pat / 100000);
	num_prf = num_age * num_zip * 2 - num_zip;				//patient ages 18 or younger do not smoke

	if (c_beta_pln != 1) {
		double coeff;
		if (c_beta_pln > 1)
			coeff = 0.8;
		else
			coeff = 1.2;
		for (int i = 0; i < beta_pln_l[0]->size(); ++i) {
			beta_pln_l[0]->at(i) *= c_beta_pln;
			beta_pln_l[1]->at(i) *= coeff;
		}
	}
	if (c_beta_qual != 1) {
		double coeff;
		if (c_beta_qual > 1)
			coeff = 0.8;
		else
			coeff = 1.2;
		for (int i = 0; i < beta_prv_nl[0]->size(); ++i) {
			beta_prv_nl[0]->at(i) *= c_beta_qual;
			beta_prv_nl[1]->at(i) *= coeff;
			beta_prv_nl[2]->at(i) *= coeff;
			beta_prv_sl[0]->at(i) *= c_beta_qual;
			beta_prv_sl[1]->at(i) *= coeff;
			beta_prv_sl[2]->at(i) *= coeff;
		}
	}
	if (c_beta_wait != 1) {
		double coeff;
		if (c_beta_wait > 1)
			coeff = 0.8;
		else
			coeff = 1.2;
		for (int i = 0; i < beta_prv_nl[0]->size(); ++i) {
			beta_prv_nl[0]->at(i) *= coeff;
			beta_prv_nl[1]->at(i) *= c_beta_wait;
			beta_prv_nl[2]->at(i) *= coeff;
			beta_prv_sl[0]->at(i) *= coeff;
			beta_prv_sl[1]->at(i) *= c_beta_wait;
			beta_prv_sl[2]->at(i) *= coeff;
		}
	}
	if (c_beta_dist != 1) {
		double coeff;
		if (c_beta_dist > 1)
			coeff = 0.8;
		else
			coeff = 1.2;
		for (int i = 0; i < beta_prv_nl[0]->size(); ++i) {
			beta_prv_nl[0]->at(i) *= coeff;
			beta_prv_nl[1]->at(i) *= coeff;
			beta_prv_nl[2]->at(i) *= c_beta_dist;
			beta_prv_sl[0]->at(i) *= coeff;
			beta_prv_sl[1]->at(i) *= coeff;
			beta_prv_sl[2]->at(i) *= c_beta_dist;
		}
	}

	cout << "*********** Creating Instance " << mods << "I" << num_pat << "_ver" << ver << " ***********" << endl;

	cout << "\nCreating Districts...\n";
	create_districts();		

	cout << "\nCreating Profiles...\n";
	create_profiles();

	cout << "\nCreating Patients...\n";
	create_patients();

	if (c_pop != 1 || c_age != 1 || c_tob != 1) {
		double slope = 0;
		int prf = 0;
		for (int d = 0; d < num_zip; ++d) {
			for (int t = 0; t < 2; ++t) {
				for (int i = 0; i < num_age; ++i) {
					
					if (i == 0 && t == 1) 
						continue;

					if (c_pop != 1) {
						if (c_pop > 1)
							slope = -0.4;
						else
							slope = 0.4;
						profiles[prf]->pop *= c_pop + slope / (num_zip - 1) * d;
					}
					if (c_tob != 1) {
						if (c_tob > 1)
							slope = 0.8;
						else
							slope = 1.2;
						if (t == 0)
							profiles[prf]->pop *= c_tob;
						else
							profiles[prf]->pop *= slope;
					}
					if (c_age != 1) {
						if (c_age > 1)
							slope = -0.4;
						else
							slope = 0.4;
						profiles[prf]->pop *= c_age + slope / (num_age - 1) * i;
					}
					++prf;
				}
			}
		}
	}


	cout << "\nCreating Providers...\n";
	create_providers();

	cout << "\nCalculating Utilities...\n";
	calc_prf_util();

	cout << "\nCreating Plans...\n";
	create_plans();
		
	cout << "\nPricing Plan 0...\n";
	price_plan0();

	cout << "\nRunning Exclusion Heuristic on Plans...\n";
	excl_heur();

	cout << "\nPricing Plans...\n";
	price_plans();

	cout << "\nWriting data to files...\n";

	write_profiles_fin();
	write_plans_fin();

	cout << "FIN\n\n\n";

	return 0;
}


void create_districts() {

	vector<double> loc_prob;
	for (int i = 0; i < num_zip; ++i) {
		vector<double> temp(2, -100);
		normal_distribution<double> zipx_dist(zip_loc[0]->at(i), zip_loc[2]->at(i));
		while (temp[0] < zip_loc[0]->at(i) - 3 * zip_loc[2]->at(i) || temp[0]> zip_loc[0]->at(i) + 3 * zip_loc[2]->at(i))
			temp[0]= round_p(zipx_dist(generator), 2);
	
		normal_distribution<double> zipy_dist(zip_loc[1]->at(i), zip_loc[3]->at(i));
		while (temp[1] < zip_loc[1]->at(i) - 3 * zip_loc[3]->at(i) || temp[1]> zip_loc[1]->at(i) + 3 * zip_loc[3]->at(i))
			temp[1] = round_p(zipy_dist(generator), 2);
		
		double min_dist = numeric_limits<double>::max();
		for (int j = 0; j < i; ++j) {
			double dist = calc_dist(temp[0] * 100 * c_dist, temp[1] * 100 * c_dist, zipx_loc[j], zipy_loc[j]);
			if (dist < min_dist)
				min_dist = dist;
		}

		if (min_dist >= 50) {
			zipx_loc.push_back(temp[0] * 100 * c_dist);
			zipy_loc.push_back(temp[1] * 100 * c_dist);
			for (int j = 0; j < i; ++j) {
				double dist = calc_dist(zipx_loc[i], zipy_loc[i], zipx_loc[j], zipy_loc[j]);
				if (dist > max_dist)
					max_dist = dist;
			}
		}
		else {
			--i;
			continue;
		}

		normal_distribution<double> pop_dist(zip_denst[0]->at(i), zip_denst[1]->at(i));
		double rand_num = -1;
		while ( rand_num < 0.1 || rand_num < zip_denst[0]->at(i) - 3 * zip_denst[1]->at(i) || rand_num > zip_denst[0]->at(i) + 3 * zip_denst[1]->at(i)) 
			rand_num = pop_dist(generator);
		
		loc_prob.push_back(round_p(rand_num, 2));

	}

	loc_dis = discrete_distribution<int>(loc_prob.begin(), loc_prob.end());	

	for (int i = 0; i < num_zip; ++i) 
		cout << fixed << setprecision(2) << "District " << i << " Zipx: " << zipx_loc[i] << " Zipy: " << zipy_loc[i] << " Pop density: " << loc_prob[i] << endl;

	cout << "\nMax Dist: " << max_dist << endl;

}

	

void create_profiles() {

	profiles.reserve(num_prf);

	for (int d = 0; d < num_zip; ++d) {
		for (int t = 0; t < 2; ++t) {
			for (int i = 0; i < num_age; ++i) {

				if (i == 0 && t == 1) 
					continue;
	
				profiles.emplace_back(new profile);
	
				profiles.back()->age = i;
				profiles.back()->zipx = zipx_loc[d];
				profiles.back()->zipy = zipy_loc[d];
				profiles.back()->tob = t;

				pat_to_prf[{i, d, t}] = profiles.size() - 1;

				for (int s = 0; s < num_ser; ++s)
					profiles.back()->pref.push_back(pat_pref[i]->at(s) * 100); 

				profiles.back()->beta_prc = -beta_pln_l[0]->at(i);
				profiles.back()->beta_brd = beta_pln_l[1]->at(i) / (util_high / 10000);
				for (int s = 0; s < num_ser; ++s) {
					if (ser_import[s] == 1) {
						profiles.back()->beta_qual.push_back(beta_prv_sl[0]->at(i));
						profiles.back()->beta_wait.push_back(-beta_prv_sl[1]->at(i));
						profiles.back()->beta_dist.push_back(-beta_prv_sl[2]->at(i));
					}
					else {
						profiles.back()->beta_qual.push_back(beta_prv_nl[0]->at(i));
						profiles.back()->beta_wait.push_back(-beta_prv_nl[1]->at(i));
						profiles.back()->beta_dist.push_back(-beta_prv_nl[2]->at(i));
					}
				}
			}
		}
	}
}


void create_patients() {

	discrete_distribution<int> age_dis(pat_age_prob.begin(), pat_age_prob.end());	

	vector<int> pop(num_zip, 0);
	for (int i = 0; i < num_pat; ++i) {

		int age = age_dis(generator);
		int distr = loc_dis(generator);

		discrete_distribution<int> tob_dis(tob_prb[age]->begin(), tob_prb[age]->end());						//0: dont use tobacco , 1: use tobacco
		int tob = tob_dis(generator);

		int prf = pat_to_prf[{age, distr, tob}];				

		++profiles[prf]->pop;
		++pop[distr];

	}

	for (int i = 0; i < num_zip; ++i) 
		cout << "District " << i << " Population: " << pop[i] << endl;

		cout << endl;

	for (int i = 0; i < num_prf; ++i) 
		cout << "Profile " << i << " population: " << profiles[i]->pop << endl;

}


void create_providers() {

	vector<int> num_prv_ser(num_ser, 0);

	uniform_real_distribution<double> ser_dist(0.0,1.0);

	vector<int> prv_type_num;
	for (int j = 0; j < prv_type_prob.size(); ++j) 
		prv_type_num.push_back(round(prv_type_prob[j] * num_prv));

	num_prv_type = prv_type_num.size();

	num_prv = 0;
	for (int j = 0; j < prv_type_prob.size(); ++j) 
		num_prv += prv_type_num[j];

	int prv_type_count = 0, prv_count = 0;

	discrete_distribution<int> qual_dis(prv_qual_prob.begin(), prv_qual_prob.end());		//Low, avg, High

	providers.reserve(num_prv);

	for (int j = 0; j < num_prv; ++j) {

		providers.emplace_back(new provider);
		providers[j]->services = vector<int>(num_ser, 0);

		if (prv_count < prv_type_num[prv_type_count]) {
			providers[j]->type = prv_type_count;
			++prv_count;
		}
		else {
			prv_count = 1;
			++prv_type_count;
			providers[j]->type = prv_type_count;
		}

		for (int s = 0; s < num_ser; ++s) {
			if (prv_ser_prob[providers[j]->type]->at(s) > 0) {
				if(prv_count <= 3 || ser_dist(generator) < prv_ser_prob[providers[j]->type]->at(s)) {
					providers[j]->services[s] = 1;
					providers[j]->num_ser += 1;
					++num_prv_ser[s];
				}
			}
		}

		tot_num_ser += exp(providers[j]->num_ser);

		int distr = loc_dis(generator);
		providers[j]->zipx = zipx_loc[distr];
		providers[j]->zipy = zipy_loc[distr];

		providers[j]->qual = qual_dis(generator) * c_qual;
	
	}

	double min_wait = 10, avg_wait = 0, temp_max = -1;

	for (int j = 0; j < num_prv; ++j) {
		for (int s = 0; s < num_ser; ++s) {
			if (providers[j]->services[s] == 1) {
				double wait = 0;
				for (int i = 0; i < num_prf; ++i)
					wait += 1.0 * profiles[i]->pop * (num_prv - num_prv_ser[s] + 1) / num_prv / (1 + calc_dist(providers[j]->zipx, providers[j]->zipy, profiles[i]->zipx, profiles[i]->zipy) / max_dist * 4);
				wait /= num_pat;
				if (wait > temp_max)
					temp_max = wait;
			}
		}
	}

	for (int j = 0; j < num_prv; ++j) {
		for (int s = 0; s < num_ser; ++s) {
			if (providers[j]->services[s] == 1) {
				double wait = 0;
				for (int i = 0; i < num_prf; ++i)
					wait += 1.0 * profiles[i]->pop * (num_prv - num_prv_ser[s] + 1) / num_prv / (1 + calc_dist(providers[j]->zipx, providers[j]->zipy, profiles[i]->zipx, profiles[i]->zipy) / max_dist * 4);
				
				wait = wait / num_pat / temp_max * 1.5;
				wait += providers[j]->qual * 0.15;

				normal_distribution<double> wait_dis(wait, 0.5);

				double rand_num = -1;	
				while (rand_num < 0 || rand_num > wait + 1)
					rand_num  = wait_dis(generator);		
				providers[j]->wait.push_back(round_p(rand_num * c_wait, 2));

				if (min_wait > providers[j]->wait[s])
					min_wait = providers[j]->wait[s];
				if (max_wait < providers[j]->wait[s])
					max_wait = providers[j]->wait[s];

				avg_wait += providers[j]->wait[s];
			}
			else 
				providers[j]->wait.push_back(0);
		}
	}

	for (int j = 0; j < num_prv; ++j) {
		for (int s = 0; s < num_ser; ++s) {
			if (providers[j]->services[s] == 1) {

				double 	rand_num = -1;
				double exp_coeff = (0.8 + 0.1 * providers[j]->wait[s] / max_wait * 4);
				double coef_s = 0.033;
				normal_distribution<double> disc_dist(exp_coeff, coef_s);
				while (rand_num < exp_coeff - 3 * coef_s || rand_num > exp_coeff + 3 * coef_s)
					rand_num  = disc_dist(generator);	
				providers[j]->neg_coeff.push_back(round_p(rand_num * c_neg, 2));
				if (same_cost)
					providers[j]->neg_coeff.back() = 1;
				rand_num = -1;

				exp_coeff = (1.8 + 0.1 * providers[j]->wait[s] / max_wait * 4);

				normal_distribution<double> self_dist(exp_coeff, coef_s);
				while (rand_num < exp_coeff - 3 * coef_s || rand_num > exp_coeff + 3 * coef_s)
					rand_num  = self_dist(generator);
				providers[j]->self_coeff.push_back(round_p(rand_num * c_sneg, 2));
			}
			else {
				providers[j]->neg_coeff.push_back(0);
				providers[j]->self_coeff.push_back(0);
			}
		}

	}

	cout << "\nMin wait: " << min_wait << " Max wait: " << max_wait << " Avg wait " << avg_wait / num_prv << endl;
	//cin.get();

	vector<double> ann_cost(num_prv, 0);

	for (int i = 0; i < num_prf; ++i) {
		double coeff = 100 / (tob_prb[profiles[i]->age]->at(0) + 1.5 * tob_prb[profiles[i]->age]->at(1));
		for (int j = 0; j < num_prv; ++j) {
			profiles[i]->pat_cost.push_back(new vector<double>(num_ser, 0));
			profiles[i]->pat_cost_self.push_back(new vector<double>(num_ser, 0));
			for (int s = 0; s < num_ser; ++s) {
				if (providers[j]->services[s] == 1) {
					profiles[i]->pat_cost[j]->at(s) = round_p(coeff * providers[j]->neg_coeff[s] * cost_ser[profiles[i]->age]->at(s) * (1 + 0.5 * profiles[i]->tob) * c_cost, 2);
					profiles[i]->pat_cost_self[j]->at(s) = round_p(coeff * providers[j]->self_coeff[s] * cost_ser[profiles[i]->age]->at(s) * (1 + 0.5 * profiles[i]->tob) * c_cost, 2);

					ann_cost[j] += profiles[i]->pat_cost[j]->at(s) * profiles[i]->pop / num_prv_ser[s] / c_neg;
				}
			}			
		}
	}


	for (int i = 0; i < num_prf; ++i) {
		for (int s = 0; s < num_ser; ++s) {
			profiles[i]->beta_wait[s] *= 4 / max_wait;
			profiles[i]->beta_dist[s] *= 4 / max_dist;
		}
	} 

	double dd = 0;
	vector<double> prv_type_ann(num_prv_type, 0);
	cout << "\nProvider types and costs:\n";
	for (int j = 0; j < num_prv; ++j) {
		providers[j]->ann_cost = 0.18 * round(ann_cost[j]) * c_anncost;
		if (providers[j]->ann_cost == 0)
			providers[j]->ann_cost = 1;
		dd += providers[j]->ann_cost;
		prv_type_ann[providers[j]->type] += providers[j]->ann_cost;
		cout << "prv " << j + 1;
		for (int pp = 0; pp <= floor(log10(num_prv) + 1) - floor(log10(j + 1) + 1); ++pp)
			cout << " ";
		cout << "type: " << providers[j]->type << " services: ";
		for (int s = 0 ; s < num_ser; ++s)
			cout << providers[j]->services[s] << " ";
		cout << " annual cost: " << providers[j]->ann_cost << endl;
	}

	if (same_ann == 1) {
		for (int j = 0; j < num_prv; ++j) 
			providers[j]->ann_cost = prv_type_ann[providers[j]->type] / prv_type_num[providers[j]->type];
	}

	cout << "TOT ANN COST PER PATIENT " << dd / num_pat << endl;

	write_profiles_raw();
	write_providers_raw();
	write_providers_fin();


}


void calc_prf_util() {

	tot_dem_prv = vector<int>(num_prv, 0);
	vector<vector<int>> dem_prv(num_prv, vector<int>(num_ser, 0));

	double min_util = numeric_limits<double>::max();
	double max_util = numeric_limits<double>::min();

	for (int i = 0; i < num_prf; ++i) {

		profiles[i]->glob_ch = vector<int>(num_ser, -1);
		profiles[i]->util_prv = vector<double>(num_prv, 0);
 
		for (int j = 0; j < num_prv; ++j) {

			profiles[i]->util_prv_ser.emplace_back(new vector<double>(num_ser, 0));
			
			double dist = calc_dist(providers[j]->zipx, providers[j]->zipy, profiles[i]->zipx, profiles[i]->zipy) / max_dist * 4;
 
			for (int s = 0; s < num_ser; ++s) {
				if (providers[j]->services[s] == 1) {
 
					profiles[i]->util_prv_ser[j]->at(s) = round(profiles[i]->beta_qual[s] * providers[j]->qual - profiles[i]->beta_dist[s] * dist - profiles[i]->beta_wait[s] * providers[j]->wait[s]);

					if (min_util > profiles[i]->util_prv_ser[j]->at(s))
						min_util = profiles[i]->util_prv_ser[j]->at(s);
					if (max_util < profiles[i]->util_prv_ser[j]->at(s))
						max_util = profiles[i]->util_prv_ser[j]->at(s);

					if (profiles[i]->glob_ch[s] == -1 || profiles[i]->util_prv_ser[j]->at(s) > profiles[i]->util_prv_ser[profiles[i]->glob_ch[s]]->at(s) || (profiles[i]->util_prv_ser[j]->at(s) == profiles[i]->util_prv_ser[profiles[i]->glob_ch[s]]->at(s) && profiles[i]->pat_cost[j]->at(s) < profiles[i]->pat_cost[profiles[i]->glob_ch[s]]->at(s))) 
						profiles[i]->glob_ch[s] = j;
				}
			}

		}


		double max_brd_util = 0;
		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser; ++s) {
				if (providers[j]->services[s] == 1) {
					if (min_util < util_low) 
						profiles[i]->util_prv_ser[j]->at(s) = profiles[i]->util_prv_ser[j]->at(s) - min_util + util_low;
					//profiles[i]->util_prv_ser[j]->at(s) = round((util_high - util_low) * (profiles[i]->util_prv_ser[j]->at(s) - min_util) / (max_util - min_util) + util_low);
 
					max_brd_util +=	profiles[i]->pref[s] * profiles[i]->util_prv_ser[j]->at(s);
				}
			}
		}

		max_brd_util = round(max_brd_util);

		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser; ++s) {
				if (providers[j]->services[s] == 1) 
					profiles[i]->util_prv[j] += round(util_high * profiles[i]->util_prv_ser[j]->at(s) * profiles[i]->pref[s] / max_brd_util);
			}
		}

		for (int s = 0; s < num_ser; ++s) {
			tot_dem_prv[profiles[i]->glob_ch[s]] += profiles[i]->pop;
			dem_prv[profiles[i]->glob_ch[s]][s] += profiles[i]->pop;
		}
	}


	cout << "\nProvider demands:\n";
	for (int j = 0; j < num_prv; ++j) {
		cout << "prv " << j + 1;
		for (int pp = 0; pp <= floor(log10(num_prv) + 1) - floor(log10(j + 1) + 1); ++pp)
			cout << " ";
		cout << "type: " << providers[j]->type << " demands:";
		for (int s = 0; s < num_ser; ++s) 
			cout << " " << dem_prv[j][s];
		cout << endl;
	}

}


void create_plans() {

	plans.emplace_back(new plan);
	for (int j = 0; j < num_prv; ++j) 								//all providers (self)
		plans.back()->providers.push_back(j);

	plans.emplace_back(new plan);
	for (int j = 0; j < num_prv; ++j) 								//all providers (insured)
		plans.back()->providers.push_back(j);

	plans.emplace_back(new plan);
	for (int j = 0; j < num_prv; ++j) {								//quality
		if (providers[j]->qual > 6) 								
			plans.back()->providers.push_back(j);
	}

	if (!plans.back()->providers.empty()) 
		plans.emplace_back(new plan);
	else
		cout << "quality plan empty\n";

	for(int j = 0; j < num_prv; ++j) {								//demand > frac_dem * num_pat
		if (tot_dem_prv[j] > 1)
			plans.back()->providers.push_back(j);
	}

	if (!plans.back()->providers.empty()) 
		plans.emplace_back(new plan);
	else
		cout << "demand plan empty\n";

	for (int j = 0; j < num_prv; ++j) {								//negotiating power (exp cost)
		double neg_coeff = 0;
		double exp_coeff = 0;
		for (int s = 0; s < num_ser; ++s) {
			if (providers[j]->services[s] == 1) {
				neg_coeff += providers[j]->neg_coeff[s];
				exp_coeff += 0.8 + 0.1 * providers[j]->wait[s] / max_wait * 4 * c_neg;
			}
		}
		if (neg_coeff <= exp_coeff) 
			plans.back()->providers.push_back(j);
	}

	if (!plans.back()->providers.empty()) 
		plans.emplace_back(new plan);
	else
		cout << "neg power plan empty\n";

	incl_heur(plans.back());									//profitibility

	if (plans.back()->providers.size() == 1) {
		delete plans.back();
		plans.pop_back();
		cout << "profitibility plan empty\n";
	}

	num_pln = plans.size();

}


void incl_heur(plan* pln) {

	pln->providers.reserve(num_prv);

	vector<int> ord_prv;
	ord_prv.reserve(num_prv);
	for (int j = 0; j < num_prv; ++j)
		ord_prv.push_back(j);
	ord_providers_cost(ord_prv, 1);

	int prv;
	for (int j = 0; j < num_prv; ++j) {
		if (providers[j]->type == 4) {
			prv = j;
			break;
		}
	}

	pln->providers.push_back(prv);
	vector<bool> prv_incl(num_prv, 0);
	prv_incl[prv] = 1;

	double cur_profit = eval_plan(pln->providers);

	bool incl_occur = 1;
	while (incl_occur) {
		incl_occur = 0;
		for (vector<int>::iterator it = ord_prv.begin(); it != ord_prv.end(); ++it) {
			if (!prv_incl[*it]) {
				pln->providers.push_back(*it);
				double incl_profit = eval_plan(pln->providers);
				if (incl_profit > cur_profit) {
					cout << "INCL RND " << " prv " << *it << " profit " << incl_profit - cur_profit << endl; 
					cur_profit = incl_profit;
					prv_incl[*it] = 1;
					incl_occur = 1;
				}
				else
					pln->providers.pop_back();
			}
		}
	}

	ord_providers(pln->providers);

}


void excl_heur() {

	for (int p = 1; p < num_pln; ++p) {
		plans.emplace_back(new plan);
		plans.back()->providers = plans[p]->providers;
	}

	for (int p = num_pln; p < plans.size(); ++p) {

		vector<bool> prv_incl(num_prv, 0);
		for (vector<int>::iterator it = plans[p]->providers.begin(); it != plans[p]->providers.end(); ++it)
			prv_incl[*it] = 1;		

		double cur_profit = eval_plan(plans[p]->providers);

		vector<int> ord_prv = plans[p]->providers;
		ord_providers_cost(ord_prv, 0);

		bool excl_occur = 1;
		int rnd = 0;
		while(excl_occur) {
			excl_occur = 0;
			for (int j = 0; j < ord_prv.size(); ++j) {
				if (ord_prv[j] != -1 && prv_incl[ord_prv[j]]) {
					int cur_prv = ord_prv[j];
					ord_prv[j] = -1;
					double excl_profit = eval_plan(ord_prv);
					if (excl_profit > cur_profit) {
						//cout << "Plan " << p << " EXCL RND " << rnd << " prv " << cur_prv << " excl prf " << excl_profit << " cur prf " << cur_profit << " profit " << excl_profit - cur_profit << endl;
						cur_profit = excl_profit;
						prv_incl[cur_prv] = 0;
						excl_occur = 1;
					}
					else
						ord_prv[j] = cur_prv;
				}
			}
			++rnd;
		}

		plans[p]->providers.clear();
		for (vector<int>::iterator it = ord_prv.begin(); it != ord_prv.end(); ++it) {
			if (*it != -1)
				plans[p]->providers.push_back(*it);
		}

		ord_providers(plans[p]->providers);

	}

	for (int p = plans.size() - 1; p >= 1; --p) {
		if (plans[p]->providers.empty()) {
			delete plans[p];
			plans.erase(plans.begin() + p);
			cout << "Red plan " << p << " empty\n";
		}
	}

	for (int p = plans.size() - 1; p >= 2; --p) {
		bool found = 0;
		for (int p2 = p - 1; p2 >= 1; --p2) {
			if (found)
				break;
			if (plans[p]->providers.size() == plans[p2]->providers.size()) {
				for (int j = 0; j < plans[p2]->providers.size(); ++j) {
					if (plans[p]->providers[j] != plans[p2]->providers[j]) 
						break;
					else if (j + 1 == plans[p]->providers.size()) {
						delete plans[p];
						plans.erase(plans.begin() + p);
						found = 1;
						cout << "Plan " << p << " duplicate of " << p2 << endl;
					}
				}
			}
		}
	}

	num_pln = plans.size();	

}


void price_plan0() {

	for (int i = 0; i < num_prf; ++i) {
		vector<int> pat_ch(num_ser, -1);
		for (vector<int>::iterator it = plans[0]->providers.begin(); it != plans[0]->providers.end(); ++it) {
			plans[0]->brd_util[i] += profiles[i]->util_prv[*it];
			for (int s = 0; s < num_ser; ++s) {
				if (providers[*it]->services[s] == 1 && (pat_ch[s] == -1 || profiles[i]->util_prv_ser[*it]->at(s) > profiles[i]->util_prv_ser[pat_ch[s]]->at(s) || (profiles[i]->util_prv_ser[*it]->at(s) == profiles[i]->util_prv_ser[pat_ch[s]]->at(s) && profiles[i]->pat_cost[*it]->at(s) < profiles[i]->pat_cost[pat_ch[s]]->at(s))))
					pat_ch[s] = *it;
			}
		}

		for (int s = 0; s < num_ser; ++s) {
			if (pat_ch[s] != -1) 
				plans[0]->cost[i] += profiles[i]->pat_cost[pat_ch[s]]->at(s);
			plans[0]->prem_comp[i] += profiles[i]->pat_cost_self[pat_ch[s]]->at(s);
		}
		
		plans[0]->cost[i] = plans[0]->prem_comp[i];
		plans[0]->util_comp[i] = profiles[i]->beta_brd * plans[0]->brd_util[i] - profiles[i]->beta_prc * plans[0]->prem_comp[i];
		plans[0]->util[i] = plans[0]->util_comp[i];
		profiles[i]->price_out = plans[0]->prem_comp[i] + 1000;
	}
}



void price_plans() {

	vector<int> max_util_pln(num_prf, 0);

	for (int i = 0; i < num_prf; ++i) {
		for (int p = 1; p < num_pln; ++p) {
			vector<int> pat_ch(num_ser, -1);
			for (vector<int>::iterator it = plans[p]->providers.begin(); it != plans[p]->providers.end(); ++it) {
				plans[p]->brd_util[i] += profiles[i]->util_prv[*it];
				for (int s = 0; s < num_ser; ++s) {
					if (providers[*it]->services[s] == 1 && (pat_ch[s] == -1 || profiles[i]->util_prv_ser[*it]->at(s) > profiles[i]->util_prv_ser[pat_ch[s]]->at(s) || (profiles[i]->util_prv_ser[*it]->at(s) == profiles[i]->util_prv_ser[pat_ch[s]]->at(s) && profiles[i]->pat_cost[*it]->at(s) < profiles[i]->pat_cost[pat_ch[s]]->at(s))))
						pat_ch[s] = *it;
				}
			}

			for (int s = 0; s < num_ser; ++s) {
				if (pat_ch[s] != -1) 
					plans[p]->cost[i] += profiles[i]->pat_cost[pat_ch[s]]->at(s);
				else 
					plans[p]->self_cost[i] += profiles[i]->pat_cost_self[profiles[i]->glob_ch[s]]->at(s);	
			}

			plans[p]->prem_comp[i] = (profiles[i]->beta_brd * (plans[p]->brd_util[i] - plans[0]->brd_util[i]) + profiles[i]->beta_prc * plans[0]->prem_comp[i]) / profiles[i]->beta_prc - plans[p]->self_cost[i] - epsilon;
			if (plans[p]->prem_comp[i] > 12000)
				plans[p]->prem_comp[i] = 12000;
			if (plans[p]->prem_comp[i] < plans[p]->cost[i]) {
				plans[p]->prem_comp[i] = profiles[i]->price_out;
				//plans[p]->prem[i] = profiles[i]->price_out;
			}
			plans[p]->util_comp[i] = profiles[i]->beta_brd * plans[p]->brd_util[i] - profiles[i]->beta_prc * (plans[p]->self_cost[i] + plans[p]->prem_comp[i]);
			if (plans[p]->util_comp[i] > plans[max_util_pln[i]]->util_comp[i])
				max_util_pln[i] = p;
		}

		for (int p = 1; p < num_pln; ++p) {
			plans[p]->prem[i] = (profiles[i]->beta_brd * plans[p]->brd_util[i] - plans[max_util_pln[i]]->util_comp[i]) / profiles[i]->beta_prc - plans[p]->self_cost[i] - 100;
			if (plans[p]->prem[i] > plans[p]->cost[i]) {
				plans[p]->single_profit += (plans[p]->prem[i] - plans[p]->cost[i]) * profiles[i]->pop;
				plans[p]->single_prem[i] = plans[p]->prem[i];
			}
			else
				plans[p]->single_prem[i] = profiles[i]->price_out;

			plans[p]->util[i] = profiles[i]->beta_brd * plans[p]->brd_util[i] - profiles[i]->beta_prc * (plans[p]->self_cost[i] + plans[p]->single_prem[i]);
		}
	}

	for (int p = 1; p < num_pln; ++p) {
		for (vector<int>::iterator it = plans[p]->providers.begin(); it != plans[p]->providers.end(); ++it) 
			plans[p]->single_profit -= providers[*it]->ann_cost;
	}


	IloEnv MP_env;
	IloModel MP_mod(MP_env);
	IloCplex MP_cpl(MP_mod);

	cout << "Building MP model\n";
	MP_cpl.setOut(MP_cpl.getEnv().getNullStream());

	build_MP(MP_cpl);

	cout << "Solving Assignment Problem\n";
	MP_cpl.solve();
	cout << "Done\n\n";

	cout << fixed << setprecision(0);

	for (int i = 0; i < num_prf; ++i) {
		for (int p = 1; p < num_pln; ++p) {
			if (MP_cpl.getValue(MP_vars.zip[i][p - 1]) < 0.5)
				plans[p]->prem[i] = profiles[i]->price_out;
			//plans[p]->util[i] = profiles[i]->beta_brd * plans[p]->brd_util[i] - profiles[i]->beta_prc * (plans[p]->self_cost[i] + plans[p]->prem[i]);
		}
	}

	double rev = 0, cost = 0, tot_rev = 0, tot_cost = 0, num_insur = 0;
	vector<bool> prv_counted_c(num_prv, 0);
	for (int i = 0; i < num_prf; ++i) {
		cout << "Prf " << i + 1;
		for (int pp = 0; pp <= floor(log10(num_prf) + 1) - floor(log10(i + 1) + 1); ++pp)
			cout << " ";
		cout << "Pop: " << profiles[i]->pop;

		for (int p = 1; p < num_pln; ++p) {
			if (MP_cpl.getValue(MP_vars.zip[i][p - 1]) > 0.5) {
				cout << "  Choice pln: " << p << "  Prem " << plans[p]->prem[i] << "  Cost: " << plans[p]->cost[i] << "  Slf cost: " << plans[p]->self_cost[i];
				double c = 0;
				for (vector<int>::iterator it = plans[p]->providers.begin(); it != plans[p]->providers.end(); ++it) 
					c += providers[*it]->ann_cost;
				cout << "  Ann pln cost: " << c;
				rev += profiles[i]->pop * plans[p]->prem[i];
				cost += profiles[i]->pop * plans[p]->cost[i];
				break;
			}
		}
		if (max_util_pln[i] > 0) {

			cout << "\nChoice pln comp: " << max_util_pln[i] << "  Prem " << plans[max_util_pln[i]]->prem_comp[i] << "  Cost: " << plans[max_util_pln[i]]->cost[i] << "  Slf cost: " << plans[max_util_pln[i]]->self_cost[i];

			tot_rev += profiles[i]->pop * plans[max_util_pln[i]]->single_prem[i];
			tot_cost += profiles[i]->pop * plans[max_util_pln[i]]->cost[i];
			for (vector<int>::iterator it = plans[max_util_pln[i]]->providers.begin(); it != plans[max_util_pln[i]]->providers.end(); ++it) {
				if (!prv_counted_c[*it]) {
					tot_cost += providers[*it]->ann_cost;
					prv_counted_c[*it] = 1;
				}
			}
			num_insur += profiles[i]->pop;
		}

		cout << endl;
	}

	cout << "\nzip:\n";
	for (int i = 0; i < num_prf; ++i) {
		cout << i + 1 << " ";
		for (int pp = 0; pp <= floor(log10(num_prf) + 1) - floor(log10(i + 1) + 1); ++pp)
			cout << " ";
		for (int p = 1; p < num_pln; ++p) 
			cout << abs(MP_cpl.getValue(MP_vars.zip[i][p - 1])) << " ";
		cout << endl;
	}

	cout << "\nxj:\n";
	for (int j = 0; j < num_prv; ++j) {
		cout << abs(MP_cpl.getValue(MP_vars.xj[j])) << " ";
		cost += providers[j]->ann_cost * MP_cpl.getValue(MP_vars.xj[j]);
	}
	cout << endl;

	double deduct = 861 * num_pat;
	
	cout << "\nNum Insured Per Mil: " << print_com(num_insur) << endl;

	cout << "\nRevenue Per Mil: " << print_com(tot_rev) << endl;
	cout << "Cost Per Mil: " << print_com(tot_cost) << endl;

	cout << "\nProfit Per Mil: " << print_com(tot_rev - tot_cost) << endl;
	cout << "\nAct Profit Per Mil: " << print_com(tot_rev - tot_cost - deduct) << endl;


	cout << "\nRevenue: " << print_com(rev) << endl;
	cout << "Cost: " << print_com(cost) << endl;

	cout << "\nProfit: " << print_com(rev - cost) << endl;
	cout << "Objective (profit): " << print_com(MP_cpl.getObjValue()) << endl;
	cout << "Actual Profit: " << print_com(rev - cost - deduct) << endl;

	

	if (MP_cpl.getStatus() != IloAlgorithm::Optimal) 
		cout << "MP NOT OPTIMAL!!!!!!!!!!!  " << MP_cpl.getStatus();

	ofstream file;
	file.open("instances.csv", std::ios::app);
	file << path_f.substr(12, path_f.length() - 13) << "," << fixed << setprecision(0) << num_pat << "," << num_prf << "," << num_prv << "," << ver << "," << deduct << "," << tot_rev << "," << tot_cost << "," << (tot_rev - tot_cost) << "," << (tot_rev - tot_cost - deduct) << "," << rev << "," << cost << "," << (rev - cost)  << "," << (rev - cost - deduct) << endl;
	file.close();



}









