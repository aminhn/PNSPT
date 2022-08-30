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

#include "utility.hpp"
#include "load_inst.hpp"
#include "build_model.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <time.h>
#include <set>
#include <iomanip>

bool plan_exists();

bool update_plans() {

	plans.emplace_back(new plan);
	plans.back()->providers.reserve(num_prv);
	for (int j = 0; j < num_prv; ++j) {
		if (SP_vars.vj_val[j] > 0.5) 
			plans.back()->providers.push_back(j);
	}

	plans.back()->prem.reserve(num_pat);
	for (int i = 0; i < num_pat; ++i) {
		if (SP_vars.zi_val[i] < 0.5)
			plans.back()->prem.push_back(patients[i]->price_out);
		else
			plans.back()->prem.push_back(SP_vars.ri_val[i]);
	}

	/*cout << "wi0s:\n";
	for (int i = 0; i < num_pat; ++i) {
		for (int s = 0; s < num_ser; ++s) 
			cout << abs(SP_vars.wi0s_val[i][s]) << " ";
		cout << endl;
	}

cin.get();*/


	++num_pln;

	calc_cost_util(plans.back());

	if (plan_exists())
		return 0;

	return 1;	

}

bool plan_exists() {

	for (int p = 0; p < num_pln - 1; ++p) {
		if (plans[num_pln - 1]->providers.size() != plans[p]->providers.size())
			continue;
		bool exists = 1;
		for (int j = 0; j < plans[num_pln - 1]->providers.size(); ++j) {
			if (plans[num_pln - 1]->providers[j] != plans[p]->providers[j]) {
				exists = 0;
				break;
			}
		}
		if (exists) {
			cout << "New plan is a copy of plan " << p << endl;
			for (int j = 0; j < plans[num_pln - 1]->providers.size(); ++j)
				cout << plans[num_pln - 1]->providers[j] << " ";
			cout << endl;
			for (int i = 0; i < num_pat; ++i)
				cout << plans[num_pln - 1]->prem[i] << " ";
			cout << endl;
			for (int i = 0; i < num_pat; ++i)
				cout << plans[p]->prem[i] << " ";

			cout << endl;
			return 1;
		}
	}

	return 0;

}


void calc_cost_util(plan* pln) {
	
	pln->brd_util.reserve(num_pat);
	pln->util.reserve(num_pat);
	pln->cost.reserve(num_pat);
	pln->self_cost.reserve(num_pat);

	//cout << setprecision(10) << fixed;

	for (int i = 0; i < num_pat; ++i) {
		
		pln->brd_util.push_back(0);
		pln->cost.push_back(0);
		pln->self_cost.push_back(0);
		pln->util.push_back(-patients[i]->beta_prc * pln->prem[i]);


		vector<int> pat_ch(num_ser, -1);
		for (vector<int>::iterator it = pln->providers.begin(); it != pln->providers.end(); ++it) {

			pln->brd_util[i] += patients[i]->util_prv[*it];
			pln->util[i] += patients[i]->beta_brd * patients[i]->util_prv[*it];

			for (int s = 0; s < num_ser; ++s) {
				if (providers[*it]->services[s] && (pat_ch[s] == -1 || patients[i]->util_prv_ser[*it]->at(s) > patients[i]->util_prv_ser[pat_ch[s]]->at(s) || (patients[i]->util_prv_ser[*it]->at(s) == patients[i]->util_prv_ser[pat_ch[s]]->at(s) && patients[i]->cost_prv[*it]->at(s) < patients[i]->cost_prv[pat_ch[s]]->at(s))))
					pat_ch[s] = *it;
			}
		}

		for (int s = 0; s < num_ser; ++s) {
			if (pat_ch[s] == -1) 
				pln->self_cost[i] += patients[i]->cost_prv_self[patients[i]->glob_ch[s]]->at(s);
			else 
				pln->cost[i] += patients[i]->cost_prv[pat_ch[s]]->at(s);
		}

		pln->util[i] -= patients[i]->beta_prc * pln->self_cost[i];

	}
}


double round_p(double num, int prec) {

	double coeff = pow(10, prec);
	return floor(num * coeff) / coeff;

}


double give_time(clock_t kk) {
	double ll = ((double)kk) / CLOCKS_PER_SEC;
	return ll;
}

string print_com(double value) {
    stringstream ss;
    ss.imbue(std::locale(""));
    ss << fixed << setprecision(0) << value;
    return ss.str();
}
