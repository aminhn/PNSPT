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
#include <numeric>      
#include <algorithm> 
#include <iomanip>


double calc_dist(double x1, double y1, double x2, double y2) {

	return round_p(sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)), 1);
}

void ord_providers_cost(vector<int>& ord_prv, bool low_to_high) {

	for (int i = 0; i < ord_prv.size(); ++i) {
		for (int j = i + 1; j < ord_prv.size(); ++j) {
			if ((low_to_high && providers[ord_prv[i]]->ann_cost > providers[ord_prv[j]]->ann_cost) || (!low_to_high && providers[ord_prv[i]]->ann_cost < providers[ord_prv[j]]->ann_cost)) {
				int temp = ord_prv[i];
				ord_prv[i] = ord_prv[j];
				ord_prv[j] = temp;
			}
		}
	}
}

void ord_providers(vector<int>& ord_prv) {

	for (int i = 0; i < ord_prv.size(); ++i) {
		for (int j = i + 1; j < ord_prv.size(); ++j) {
			if (ord_prv[i] > ord_prv[j]) {
				int temp = ord_prv[i];
				ord_prv[i] = ord_prv[j];
				ord_prv[j] = temp;
			}
		}
	}
}

double round_p(double num, int prec) {

	double coeff = pow(10, prec);

	return floor(num * coeff) / coeff;

}

double eval_plan(vector<int>& prvs) {

	double profit = 0;
	for (vector<int>::iterator it = prvs.begin(); it != prvs.end(); ++it) {
		if (*it != -1)
			profit -= providers[*it]->ann_cost;
	}

	for (int i = 0; i < num_prf; ++i) {
		vector<int> pat_ch(num_ser, -1);
		double net_util = 0, cost = 0, self_cost = 0;
		for (vector<int>::iterator it = prvs.begin(); it != prvs.end(); ++it) {
			if (*it == -1)
				continue;
			net_util += profiles[i]->util_prv[*it];
			for (int s = 0; s < num_ser; ++s) {
				if (providers[*it]->services[s] == 1 && (pat_ch[s] == -1 || profiles[i]->util_prv_ser[*it]->at(s) > profiles[i]->util_prv_ser[pat_ch[s]]->at(s) || (profiles[i]->util_prv_ser[*it]->at(s) == profiles[i]->util_prv_ser[pat_ch[s]]->at(s) && profiles[i]->pat_cost[*it]->at(s) < profiles[i]->pat_cost[pat_ch[s]]->at(s))))
					pat_ch[s] = *it;
			}
		}

		for (int s = 0; s < num_ser; ++s) {
			if (pat_ch[s] != -1) 
				cost += profiles[i]->pat_cost[pat_ch[s]]->at(s);
			else 
				self_cost += profiles[i]->pat_cost_self[profiles[i]->glob_ch[s]]->at(s);	
		}

		double prem = (profiles[i]->beta_brd * (net_util - plans[0]->brd_util[i]) + profiles[i]->beta_prc * plans[0]->prem_comp[i]) / profiles[i]->beta_prc - self_cost - epsilon;
		if (prem > 12000)
			prem = 12000;
		if (prem > cost) 
			profit += (prem - cost) * profiles[i]->pop;
	}

	return profit;

}

string print_com(double value) {
    stringstream ss;
    ss.imbue(std::locale(""));
    ss << fixed << setprecision(0) << value;
    return ss.str();
}






