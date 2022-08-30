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

#include "load_inst.hpp"
#include "utility.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <string>

int num_pat, num_prv, num_pln, num_plnc, num_ser, tot_num_pat = 0, max_single_pln, epsilon = 1900;
double max_single_prof = 0, deduct;

vector<patient*> patients;
vector<provider*> providers;
vector<plan*> plans;
vector<plan*> comp_plans;

bool load_patients();
bool load_providers();
bool load_plans();

bool load_inst() {

	if (!load_patients())
		return 0;
	num_pat = patients.size();

	if (!load_providers())
		return 0;
	num_prv = providers.size();
	num_ser = providers.back()->services.size();

	if (!load_plans())
		return 0;
	num_pln = plans.size();
	num_plnc = comp_plans.size();

	for (int i = 0; i < num_pat; ++i)
		tot_num_pat += patients[i]->pop;

	deduct = 861 * tot_num_pat;

	cout << "\n\nnum_pat: " << num_pat << ", num_prv: " << num_prv << ", num_ser: " << num_ser << ", num_pln: " << num_pln << ", num_plnc: " << num_plnc << endl << endl;
	cout << "All data loaded\n";

	return 1;
}


bool load_patients() {

	string in_file = in_path + "prf.csv";

	cout << "Loading " << in_file << endl;
	ifstream file(in_file);
	if (file.good()) {
		string line;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			patients.emplace_back(new patient);
			istringstream word(line);
			string itm;
			int pos = 0;
			while (getline(word, itm, ',')) {
				if (pos == 0)
					patients.back()->pop = stoi(itm);
				else 
					patients.back()->price_out = stoi(itm);
				++pos;
			}
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();

	in_file = in_path + "beta_pln.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			int pos = 0;
			while (getline(word, itm, ',')) {
				if (pos == 0)
					patients[i]->beta_brd = abs(stod(itm));
				else
					patients[i]->beta_prc = abs(stod(itm));
				++pos;
			}
			++i;
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();

	
	in_file = in_path + "glob_ch.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			while (getline(word, itm, ',')) 
				patients[i]->glob_ch.push_back(stoi(itm));
			++i;
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();

	in_file = in_path + "prf_pref.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			while (getline(word, itm, ',')) 
				patients[i]->pref.push_back(stod(itm));
			++i;
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();



	in_file = in_path + "prf_util_prv.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			while (getline(word, itm, ',')) 
				patients[i]->util_prv.push_back(stod(itm));
			++i;
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();


	in_file = in_path + "prf_util_prv_ser.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			patients[i]->util_prv_ser.push_back(new vector<double>);
			istringstream word(line);
			string itm;
			if (!isdigit(line[0]) && line.size() == 2) {
				++i;
				continue;
			}
			while (getline(word, itm, ',')) 
				patients[i]->util_prv_ser.back()->push_back(stod(itm));
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();

	in_file = in_path + "prf_cost.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			patients[i]->cost_prv.push_back(new vector<double>);
			istringstream word(line);
			string itm;
			if (!isdigit(line[0]) && line.size() == 2) {
				++i;
				continue;
			}
			while (getline(word, itm, ',')) 
				patients[i]->cost_prv.back()->push_back(stod(itm));
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();

	in_file = in_path + "prf_cost_self.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			patients[i]->cost_prv_self.push_back(new vector<double>);
			istringstream word(line);
			string itm;
			if (!isdigit(line[0]) && line.size() == 2) {
				++i;
				continue;
			}
			while (getline(word, itm, ',')) 
				patients[i]->cost_prv_self.back()->push_back(stod(itm));
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();


	return 1;
}


bool load_providers() {

	string in_file = in_path + "prv.csv";

	cout << "Loading " << in_file << endl;
	ifstream file(in_file);
	if (file.good()) {
		string line;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			providers.emplace_back(new provider);
			istringstream word(line);
			string itm;
			int pos = 0;
			providers.back()->services.reserve(num_ser);
			providers.back()->num_ser = 0;
			while (getline(word, itm, ',')) {
				if (pos == 0)
					providers.back()->ann_cost = stol(itm);
				else {
					providers.back()->services.push_back(stoi(itm) == 1);
					providers.back()->num_ser += stoi(itm);
				}
				++pos;
			}
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();

	return 1;
}


bool load_plans() {

	string in_file = in_path + "pln.csv";

	cout << "Loading " << in_file << endl;
	ifstream file(in_file);
	if (file.good()) {
		string line;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			comp_plans.emplace_back(new plan);
			istringstream word(line);
			string itm;
			while (getline(word, itm, ','))
				comp_plans.back()->providers.push_back(stoi(itm));
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();



	for (int p = 1; p < comp_plans.size(); ++p) {
		plans.push_back(new plan);
		plans.back()->providers = comp_plans[p]->providers;	
	}


	in_file = in_path + "epsilon.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			while (getline(word, itm, ',')) 
				epsilon = stoi(itm);
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
	}
	file.close();


	in_file = in_path + "single_pln_prof.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int p = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			while (getline(word, itm, ',')) {
				if (stod(itm) > max_single_prof) {
					max_single_prof = stod(itm);
					max_single_pln = p;
				}
			}
			++p;
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();




	in_file = in_path + "pln_prem_comp.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			comp_plans[i]->prem.reserve(num_pat);
			while (getline(word, itm, ',')) 
				comp_plans[i]->prem.push_back(stod(itm));
			++i;
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();

	in_file = in_path + "pln_prem.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			plans[i]->prem.reserve(num_pat);
			while (getline(word, itm, ',')) 
				plans[i]->prem.push_back(stod(itm));
			++i;
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();

	in_file = in_path + "pln_prem_single.csv";
	cout << "Loading " << in_file << endl;
	file.open(in_file);
	if (file.good()) {
		string line;
		int i = 0;
		while (getline(file, line)) {
			if (!isdigit(line[0]) && (line.size() == 1 || !isdigit(line[1])))
				continue;
			istringstream word(line);
			string itm;
			plans[i]->single_prem.reserve(num_pat);
			while (getline(word, itm, ',')) 
				plans[i]->single_prem.push_back(stod(itm));
			++i;
		}
	}
	else {
		cout << "!!!!!! No such file exists: " << in_file << " !!!!!!\n";
		return 0;
	}
	file.close();


	//for (int p = plans.size() - 1; p > 0; --p)
	//	plans.pop_back();


	for (int p = 0; p < plans.size(); ++p) 
		calc_cost_util(plans[p]);

	for (int p = 0; p < comp_plans.size(); ++p)
		calc_cost_util(comp_plans[p]);
	


	return 1;
}

