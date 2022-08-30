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


void write_profiles_raw() {

	string out_file = path_r + "prf.csv";

	ofstream file;
	file.open(out_file);
	file << "prf,pop,age,zipx,zipy,tob\n";
	for (int i = 0; i < num_prf; ++i) {
		file << i + 1 << "," << profiles[i]->pop << "," << profiles[i]->age << "," << profiles[i]->zipx << "," << profiles[i]->zipy << "," << profiles[i]->tob << endl;
	}
	file.close();

}

void write_profiles_fin() {

	string out_file = path_f + "prf.csv";
	ofstream file;
	file.open(out_file);
	file << "population,priceout\n";
	for (int i = 0; i < num_prf; ++i) 
		file << profiles[i]->pop << "," << profiles[i]->price_out << endl;
	file.close();

	out_file = path_f + "glob_ch.csv";
	file.open(out_file);
	for (int i = 0; i < num_prf; ++i) {
		for (int s = 0; s < num_ser - 1; ++s) 
			file << profiles[i]->glob_ch[s] << ",";
		file << profiles[i]->glob_ch.back() << endl;
	}
	file.close();

	out_file = path_f + "prf_cost.csv";
	file.open(out_file);
	for (int i = 0; i < num_prf; ++i) {
		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser - 1; ++s)
				file << profiles[i]->pat_cost[j]->at(s) << ",";
			file << profiles[i]->pat_cost[j]->back() << endl;
		}
		file << "-1\n";
	}
	file.close();

	out_file = path_f + "prf_cost_self.csv";
	file.open(out_file);
	for (int i = 0; i < num_prf; ++i) {
		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser - 1; ++s)
				file << profiles[i]->pat_cost_self[j]->at(s) << ",";
			file << profiles[i]->pat_cost_self[j]->back() << endl;
		}
		file << "-1\n";
	}
	file.close();

	out_file = path_f + "beta_pln.csv";
	file.open(out_file);
	for (int i = 0; i < num_prf; ++i) 
		file << profiles[i]->beta_brd << "," << profiles[i]->beta_prc << endl;
	file.close();

	out_file = path_f + "prf_util_prv_ser.csv";
	file.open(out_file);
	for (int i = 0; i < num_prf; ++i) {
		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser - 1; ++s)
				file << profiles[i]->util_prv_ser[j]->at(s) << ",";
			file << profiles[i]->util_prv_ser[j]->back() << endl;
		}
		file << "-1\n";
	}
	file.close();


	out_file = path_f + "prf_util_prv.csv";
	file.open(out_file);
	for (int i = 0; i < num_prf; ++i) {
		for (int j = 0; j < num_prv - 1; ++j)
			file << profiles[i]->util_prv[j] << ",";
		file << profiles[i]->util_prv.back() << endl;
	}
	file.close();

	out_file = path_f + "prf_pref.csv";
	file.open(out_file);
	for (int i = 0; i < num_prf; ++i) {
		for (int s = 0; s < num_ser - 1; ++s) 
			file << fixed << setprecision(7) << profiles[i]->pref[s] << ",";
		file << fixed << setprecision(7) << profiles[i]->pref.back() << endl;
	}
	file.close();

}


void write_providers_raw() {

	string out_file = path_r + "prv.csv";

	ofstream file;
	file.open(out_file);

	file << "type,zipx,zipy,quality,annual_cost,wait,negotiation_coeff,self_coeff,services" << endl;

	for (int j = 0; j < providers.size(); ++j) {
		file << providers[j]->type << "," << providers[j]->zipx << "," << providers[j]->zipy << "," << providers[j]->qual << "," << providers[j]->ann_cost;
		for (int s = 0; s < num_ser; ++s)
			file << "," << providers[j]->wait[s];
		for (int s = 0; s < num_ser; ++s)
			file << "," << providers[j]->neg_coeff[s];
		for (int s = 0; s < num_ser; ++s)
			file << "," << providers[j]->self_coeff[s];
		for (int s = 0; s < num_ser; ++s)
			file << "," << providers[j]->services[s];

		file << endl;
	}
	file.close();
}

void write_providers_fin() {

	string out_file = path_f + "prv.csv";

	ofstream file;
	file.open(out_file);

	file << "annual_cost,services" << endl;

	for (int j = 0; j < providers.size(); ++j) {
		file << fixed << setprecision(0) << providers[j]->ann_cost;
		for (int s = 0; s < num_ser; ++s)
			file << "," << providers[j]->services[s];
		file << endl;
	}
	file.close();
}

void write_plans_fin() {

	string out_file = path_f + "pln.csv";
	ofstream file;
	file.open(out_file);
	for (int p = 0; p < plans.size(); ++p) {
		for (int j = 0; j < plans[p]->providers.size() - 1; ++j)
			file << plans[p]->providers[j] << ",";
		file << plans[p]->providers.back() << endl;
	}
	file.close();

	out_file = path_f + "epsilon.csv";
	file.open(out_file);
	file << epsilon << endl;
	file.close();

	out_file = path_f + "pln_prem.csv";
	file.open(out_file);
	for (int p = 1; p < num_pln; ++p) {
		for (int i = 0; i < num_prf - 1; ++i) 
			file << plans[p]->prem[i] << ",";
		file << plans[p]->prem.back() << endl;
	}
	file.close();

	out_file = path_f + "pln_prem_comp.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_prf - 1; ++i) 
			file << plans[p]->prem_comp[i] << ",";
		file << plans[p]->prem_comp.back() << endl;
	}
	file.close();

	out_file = path_f + "pln_prem_single.csv";
	file.open(out_file);
	for (int p = 1; p < num_pln; ++p) {
		for (int i = 0; i < num_prf - 1; ++i) 
			file << plans[p]->single_prem[i] << ",";
		file << plans[p]->single_prem.back() << endl;
	}
	file.close();


	out_file = path_f + "pln_cost.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_prf - 1; ++i) 
			file << plans[p]->cost[i] << ",";
		file << plans[p]->cost.back() << endl;
	}
	file.close();

	out_file = path_f + "pln_self_cost.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_prf - 1; ++i) 
			file << plans[p]->self_cost[i] << ",";
		file << plans[p]->self_cost.back() << endl;
	}
	file.close();

	out_file = path_f + "pln_brd_util.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_prf - 1; ++i) 
			file << fixed << setprecision(0) << plans[p]->brd_util[i] << ",";
		file << plans[p]->brd_util.back() << endl;
	}
	file.close();

	out_file = path_f + "pln_util.csv";
	file.open(out_file);
	for (int p = 1; p < num_pln; ++p) {
		for (int i = 0; i < num_prf - 1; ++i) 
			file << fixed << setprecision(0) << plans[p]->util[i] << ",";
		file << plans[p]->util.back() << endl;
	}
	file.close();

	out_file = path_f + "pln_util_comp.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_prf - 1; ++i) 
			file << fixed << setprecision(0) << plans[p]->util_comp[i] << ",";
		file << plans[p]->util_comp.back() << endl;
	}
	file.close();

	out_file = path_f + "single_pln_prof.csv";
	file.open(out_file);
	for (int p = 1; p < num_pln; ++p) {
		file << fixed << setprecision(0) << plans[p]->single_profit << endl;
	}
	file.close();

} 

bool load_data(string inst, vector<vector<int>*>& dat) {

	ifstream file(inst);
	if (file.good()) {
		string line;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			dat.push_back(new vector<int>);
			while (getline(word, itm, ','))
				dat.back()->push_back(stoi(itm));
		}

	}
	else {
		if (!file.good())
			cout << "!!!!!! No such file exists: " << inst << " !!!!!!\n";
		return 0;
	}

	return 1;
}


bool load_data(string inst, vector<vector<double>*>& dat) {

	//cout << "Loading " << inst << endl;
	ifstream file(inst);
	if (file.good()) {
		string line;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			dat.push_back(new vector<double>);
			while (getline(word, itm, ','))
				dat.back()->push_back(stod(itm));
		}

	}
	else {
		if (!file.good())
			cout << "!!!!!! No such file exists: " << inst << " !!!!!!\n";
		return 0;
	}

	//cout << "DONE\n";

	return 1;
}

bool load_data(string inst, vector<double>& dat) {

	//cout << "Loading " << inst << endl;
	ifstream file(inst);
	if (file.good()) {
		string line;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			while (getline(word, itm, ','))
				dat.push_back(stod(itm));
		}
	}
	else {
		if (!file.good())
			cout << "!!!!!! No such file exists: " << inst << " !!!!!!\n";
		return 0;
	}

	//cout << "DONE\n";

	return 1;
}

bool load_data(string inst, vector<int>& dat) {

	//cout << "Loading " << inst << endl;
	ifstream file(inst);
	if (file.good()) {
		string line;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			while (getline(word, itm, ','))
				dat.push_back(stoi(itm));
		}
	}
	else {
		if (!file.good())
			cout << "!!!!!! No such file exists: " << inst << " !!!!!!\n";
		return 0;
	}

	//cout << "DONE\n";

	return 1;
}
