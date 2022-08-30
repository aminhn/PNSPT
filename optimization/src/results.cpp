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

#include "build_model.hpp"
#include "utility.hpp"
#include <set>

vector<int> num_self_insur(3, 0); 		//pos 0: only one plan offered by pay, pos 1: multiple plans offered by pay, pos 2: optimising mulitple plans (post opt) by pay
vector<int> num_comp_insur(3, 0);
vector<int> num_pay_insur(3, 0);
vector<int> num_tot_insur(3, 0);
vector<int> num_pln_off(3, 0);
vector<int> num_pln_off_comp(3, 0);
vector<int> pay_prv_cont(3, 0);
vector<int> comp_prv_cont(3, 0);
vector<int> prv_tot_cont(3, 0);
vector<double> pat_tot_pay(3, 0);
vector<double> pat_tot_cost(3, 0);
vector<double> prv_fix_cost(3, 0);
vector<double> pay_prof(3, 0);
vector<double> comp_prof(3, 0);
vector<double> pay_act_prof(3, 0);
vector<double> comp_act_prof(3, 0);
vector<double> pay_prem_rec(3, 0);
vector<double> comp_prem_rec(3, 0);
vector<double> pay_pat_cost(3, 0);
vector<double> comp_pat_cost(3, 0);
vector<double> self_pat_cost(3, 0);
vector<double> pay_fix_cost(3, 0);
vector<double> comp_fix_cost(3, 0);
vector<double> avg_brd_util(3, 0);
vector<double> avg_util(3, 0);

string out_path;

void write_plans();
void write_model_out(IloCplex& IP_cpl, bool fin);
void write_stats();
void rec_single_pln();
double diff(vector<int>& vec, int pos1, int denom);
double diff(vector<int>& vec, int pos1);
double diff(vector<double>& vec, int pos1, int denom);
double diff(vector<double>& vec, int pos1);

void rec_results(IloCplex& IP_cpl, bool fin) {

	out_path = "./results/" + instance + "/";

	if(!fin) {
		string sys_exc1= "rm -rf " + out_path;
		string sys_exc2= "mkdir -p " + out_path;
		bool exc = system(sys_exc1.c_str());
		exc = system(sys_exc2.c_str());
	}

	write_model_out(IP_cpl, fin);

	int pos = 1;

	if (fin) {
		write_plans();
		pos = 2;
	}
	else
		rec_single_pln();

			
	vector<int> pln_off(num_pln, 0);
	vector<int> pln_off_comp(num_plnc, 0);
	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_pln; ++p) {
			if (IP_cpl.getValue(MP_vars.yip[i][p]) > 0.5) {
				pln_off[p] = 1;
				num_pay_insur[pos] += patients[i]->pop;
				num_tot_insur[pos] += patients[i]->pop;
				pat_tot_pay[pos] += (plans[p]->prem[i] + plans[p]->self_cost[i]) * patients[i]->pop;
				pat_tot_cost[pos] += (plans[p]->cost[i] + plans[p]->self_cost[i]) * patients[i]->pop;
				pay_prof[pos] += (plans[p]->prem[i] - plans[p]->cost[i]) * patients[i]->pop;
				pay_prem_rec[pos] += plans[p]->prem[i] * patients[i]->pop;
				pay_pat_cost[pos] += plans[p]->cost[i] * patients[i]->pop;
				self_pat_cost[pos] += plans[p]->self_cost[i] * patients[i]->pop;
				avg_brd_util[pos] += plans[p]->brd_util[i] * patients[i]->pop;
				avg_util[pos] += plans[p]->util[i] * patients[i]->pop;
				break;
			}	
		}

		for (int p = 0; p < num_plnc; ++p) {
			if (IP_cpl.getValue(MP_vars.yip_comp[i][p]) > 0.5) {
				if (p == 0) {
					num_self_insur[pos] += patients[i]->pop;
					pat_tot_pay[pos] += comp_plans[p]->cost[i] * patients[i]->pop;
					pat_tot_cost[pos] += comp_plans[p]->cost[i] * patients[i]->pop;
					self_pat_cost[pos] += comp_plans[p]->cost[i] * patients[i]->pop;
					avg_brd_util[pos] += comp_plans[p]->brd_util[i] * patients[i]->pop;
					avg_util[pos] += comp_plans[p]->util[i] * patients[i]->pop;
					break;
				}
				else {
					pln_off_comp[p] = 1;
					num_comp_insur[pos] += patients[i]->pop;
					num_tot_insur[pos] += patients[i]->pop;
					pat_tot_pay[pos] += (comp_plans[p]->prem[i] + comp_plans[p]->self_cost[i]) * patients[i]->pop;
					pat_tot_cost[pos] += (comp_plans[p]->cost[i] + comp_plans[p]->self_cost[i]) * patients[i]->pop;
					comp_prof[pos] += (comp_plans[p]->prem[i] - comp_plans[p]->cost[i]) * patients[i]->pop;
					comp_prem_rec[pos] += comp_plans[p]->prem[i] * patients[i]->pop;
					comp_pat_cost[pos] += comp_plans[p]->cost[i] * patients[i]->pop;
					self_pat_cost[pos] += comp_plans[p]->self_cost[i] * patients[i]->pop;
					avg_brd_util[pos] += comp_plans[p]->brd_util[i] * patients[i]->pop;
					avg_util[pos] += comp_plans[p]->util[i] * patients[i]->pop;
					break;
				}
			}
		}
	}
	
	for (int j = 0; j < num_prv; ++j) {
		if (IP_cpl.getValue(MP_vars.xj[j]) > 0.5) {
			++prv_tot_cont[pos];
			++pay_prv_cont[pos];
			pay_fix_cost[pos] += providers[j]->ann_cost;
			prv_fix_cost[pos] += providers[j]->ann_cost;
			pay_prof[pos] -= providers[j]->ann_cost;
		}
	}

	pay_act_prof[pos] = pay_prof[pos] - deduct;
	avg_brd_util[pos] /= tot_num_pat;
	avg_util[pos] /= tot_num_pat;

	for (int p = 0; p < num_pln; ++p) 
		num_pln_off[pos] += pln_off[p];

	vector<int> prv_cont_comp(num_prv, 0);
	for (int p = 1; p < num_plnc; ++p) {
		if (pln_off_comp[p] == 1) {
			++num_pln_off_comp[pos];
			for (vector<int>::iterator it = comp_plans[p]->providers.begin(); it != comp_plans[p]->providers.end(); ++it)
				prv_cont_comp[*it] = 1;
		}
	}

	for (int j = 0; j < num_prv; ++j) {
		if (prv_cont_comp[j] == 1) {
			++prv_tot_cont[pos];
			++comp_prv_cont[pos];
			comp_fix_cost[pos] += providers[j]->ann_cost;
			prv_fix_cost[pos] += providers[j]->ann_cost;
			comp_prof[pos] -= providers[j]->ann_cost;
		}
	} 

	comp_act_prof[pos] = comp_prof[pos] - deduct;

	if (fin) 
		write_stats();
		
}


void rec_single_pln() {


	vector<int> max_util_pln(num_pat, 0);
	for (int i = 0; i < num_pat; ++i) {
		for (int p = 1; p < num_plnc; ++p) {
			if (comp_plans[p]->util[i] > comp_plans[max_util_pln[i]]->util[i]) 
				max_util_pln[i] = p;
		}
	}

	vector<int> pln_off_comp(num_plnc, 0);
	for (int i = 0; i < num_pat; ++i) {
		if (plans[max_single_pln]->single_prem[i] != patients[i]->price_out) {
			num_pay_insur[0] += patients[i]->pop;
			num_tot_insur[0] += patients[i]->pop;
			pat_tot_pay[0] += (plans[max_single_pln]->single_prem[i] + plans[max_single_pln]->self_cost[i]) * patients[i]->pop;
			pat_tot_cost[0] += (plans[max_single_pln]->cost[i] + plans[max_single_pln]->self_cost[i]) * patients[i]->pop;
			pay_prof[0] += (plans[max_single_pln]->single_prem[i] - plans[max_single_pln]->cost[i]) * patients[i]->pop;
			pay_prem_rec[0] += plans[max_single_pln]->single_prem[i] * patients[i]->pop;
			pay_pat_cost[0] += plans[max_single_pln]->cost[i] * patients[i]->pop;
			self_pat_cost[0] += plans[max_single_pln]->self_cost[i] * patients[i]->pop;
			avg_brd_util[0] += plans[max_single_pln]->brd_util[i] * patients[i]->pop;
			avg_util[0] += plans[max_single_pln]->util[i] * patients[i]->pop;
		}	
		else if (max_util_pln[i] == 0) {
			num_self_insur[0] += patients[i]->pop;
			pat_tot_pay[0] += comp_plans[0]->cost[i] * patients[i]->pop;
			pat_tot_cost[0] += comp_plans[0]->cost[i] * patients[i]->pop;
			self_pat_cost[0] += comp_plans[0]->cost[i] * patients[i]->pop;
			avg_brd_util[0] += comp_plans[0]->brd_util[i] * patients[i]->pop;
			avg_util[0] += comp_plans[0]->util[i] * patients[i]->pop;
		}
		else {
			pln_off_comp[max_util_pln[i]] = 1;
			num_comp_insur[0] += patients[i]->pop;
			num_tot_insur[0] += patients[i]->pop;
			pat_tot_pay[0] += (comp_plans[max_util_pln[i]]->prem[i] + comp_plans[max_util_pln[i]]->self_cost[i]) * patients[i]->pop;
			pat_tot_cost[0] += (comp_plans[max_util_pln[i]]->cost[i] + comp_plans[max_util_pln[i]]->self_cost[i]) * patients[i]->pop;
			comp_prof[0] += (comp_plans[max_util_pln[i]]->prem[i] - comp_plans[max_util_pln[i]]->cost[i]) * patients[i]->pop;
			comp_prem_rec[0] += comp_plans[max_util_pln[i]]->prem[i] * patients[i]->pop;
			comp_pat_cost[0] += comp_plans[max_util_pln[i]]->cost[i] * patients[i]->pop;
			self_pat_cost[0] += comp_plans[max_util_pln[i]]->self_cost[i] * patients[i]->pop;
			avg_brd_util[0] += comp_plans[max_util_pln[i]]->brd_util[i] * patients[i]->pop;
			avg_util[0] += comp_plans[max_util_pln[i]]->util[i] * patients[i]->pop;
		}
	}

	for (vector<int>::iterator it = plans[max_single_pln]->providers.begin(); it != plans[max_single_pln]->providers.end(); ++it) {
		++prv_tot_cont[0];
		++pay_prv_cont[0];
		pay_fix_cost[0] += providers[*it]->ann_cost;
		prv_fix_cost[0] += providers[*it]->ann_cost;
		pay_prof[0] -= providers[*it]->ann_cost;	
	}

	pay_act_prof[0] = pay_prof[0] - deduct;
	avg_brd_util[0] /= tot_num_pat;
	avg_util[0] /= tot_num_pat;

	num_pln_off[0] = 1;

	vector<int> prv_cont_comp(num_prv, 0);
	for (int p = 1; p < num_plnc; ++p) {
		if (pln_off_comp[p] == 1) {
			++num_pln_off_comp[0];
			for (vector<int>::iterator it = comp_plans[p]->providers.begin(); it != comp_plans[p]->providers.end(); ++it)
				prv_cont_comp[*it] = 1;
		}
	}

	for (int j = 0; j < num_prv; ++j) {
		if (prv_cont_comp[j] == 1) {
			++prv_tot_cont[0];
			++comp_prv_cont[0];
			comp_fix_cost[0] += providers[j]->ann_cost;
			prv_fix_cost[0] += providers[j]->ann_cost;
			comp_prof[0] -= providers[j]->ann_cost;
		}
	}

	comp_act_prof[0] = comp_prof[0] - deduct;

}


void write_model_out(IloCplex& IP_cpl, bool fin) {

	string out_path2;

	if (fin)
		out_path2 = out_path + "fin_";
	else
		out_path2 = out_path + "org_";

	string out_file = out_path2 + "zip.csv";
	ofstream file;
	file.open(out_file);
	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_pln - 1; ++p) 
			file << IP_cpl.getValue(MP_vars.zip[i][p]) << ",";
		file << IP_cpl.getValue(MP_vars.zip[i][num_pln - 1]) << endl;
	}
	file.close();

	out_file = out_path2 + "yip_comp.csv";
	file.open(out_file);
	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_plnc - 1; ++p) 
			file << IP_cpl.getValue(MP_vars.yip_comp[i][p]) << ",";
		file << IP_cpl.getValue(MP_vars.yip_comp[i][num_plnc - 1]) << endl;
	}
	file.close();

	out_file = out_path2 + "yip.csv";
	file.open(out_file);
	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_pln - 1; ++p) 
			file << IP_cpl.getValue(MP_vars.yip[i][p]) << ",";
		file << IP_cpl.getValue(MP_vars.yip[i][num_pln - 1]) << endl;
	}
	file.close();

	out_file = out_path2 + "xj.csv";
	file.open(out_file);
	for (int j = 0; j < num_prv - 1; ++j) 
		file << IP_cpl.getValue(MP_vars.xj[j]) << ",";
	file << IP_cpl.getValue(MP_vars.xj[num_prv - 1]) << endl;
	file.close();

}


void write_plans() {

	string out_file = out_path + "pln.csv";
	ofstream file;
	file.open(out_file);
	for (int p = 0; p < plans.size(); ++p) {
		for (int j = 0; j < plans[p]->providers.size() - 1; ++j)
			file << plans[p]->providers[j] << ",";
		file << plans[p]->providers.back() << endl;
	}
	file.close();

	out_file = out_path + "pln_prem.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_pat - 1; ++i) 
			file << plans[p]->prem[i] << ",";
		file << plans[p]->prem.back() << endl;
	}
	file.close();

	out_file = out_path + "pln_cost.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_pat - 1; ++i) 
			file << plans[p]->cost[i] << ",";
		file << plans[p]->cost.back() << endl;
	}
	file.close();

	out_file = out_path + "pln_self_cost.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_pat - 1; ++i) 
			file << plans[p]->self_cost[i] << ",";
		file << plans[p]->self_cost.back() << endl;
	}
	file.close();

	out_file = out_path + "pln_brd_util.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_pat - 1; ++i) 
			file << plans[p]->brd_util[i] << ",";
		file << plans[p]->brd_util.back() << endl;
	}
	file.close();


	out_file = out_path + "pln_util.csv";
	file.open(out_file);
	for (int p = 0; p < num_pln; ++p) {
		for (int i = 0; i < num_pat - 1; ++i) 
			file << plans[p]->util[i] << ",";
		file << plans[p]->util.back() << endl;
	}
	file.close();

}


void write_stats () {

	string res_file = "./results.csv";

	ofstream file;
	file.open(res_file, std::ios::app);
	for (int f = 0; f < 2; ++f) {
		file << instance << ",";
		if (f == 0)
			file << "single vs multiple ";
		else
			file << "pre-opt vs post-opt ";

		file << fixed << setprecision(2) << "percentages:," << diff(num_self_insur, f, 0) << "," << diff(num_comp_insur, f, 0) << "," << diff(num_pay_insur, f, 0) << "," << diff(num_tot_insur, f, 0) 
		<< "," << diff(num_pln_off_comp, f, 0) << "," << diff(num_pln_off, f, 0) << "," << diff(pat_tot_pay, f, 0) << "," << diff(pat_tot_cost, f, 0) << "," << diff(self_pat_cost, f, 0) << "," 
		<< diff(comp_pat_cost, f, 0) << "," << diff(pay_pat_cost, f, 0) << "," << diff(prv_tot_cont, f, num_prv) << "," << diff(comp_prv_cont, f, num_prv) << "," << diff(pay_prv_cont, f, num_prv) << "," << diff(prv_fix_cost, f, 0)
 		<< "," << diff(comp_fix_cost, f, 0) << "," << diff(pay_fix_cost, f, 0) << "," << diff(comp_prem_rec, f, 0) << "," << diff(pay_prem_rec, f, 0) << "," << diff(comp_prof, f, 0) << "," << diff(pay_prof, f, 0)  << "," 
		<< diff(comp_act_prof, f, 0) << "," << diff(pay_act_prof, f, 0) << "," << diff(avg_brd_util, f, 0) << "," << diff(avg_util, f, 0) << endl;
	
		file << instance << ",";
		if (f == 0)
			file << "single vs multiple ";
		else
			file << "pre-opt vs post-opt ";

		file << fixed << setprecision(0) << "actuals:," << diff(num_self_insur, f) << "," << diff(num_comp_insur, f) << "," << diff(num_pay_insur, f) << "," << diff(num_tot_insur, f) 
		<< "," << diff(num_pln_off_comp, f) << "," << diff(num_pln_off, f) << "," << diff(pat_tot_pay, f) << "," << diff(pat_tot_cost, f) << "," << diff(self_pat_cost, f) << "," << diff(comp_pat_cost, f) << "," << diff(pay_pat_cost, f)
		<< "," << diff(prv_tot_cont, f) << "," << diff(comp_prv_cont, f) << "," << diff(pay_prv_cont, f) << "," << diff(prv_fix_cost, f) << "," << diff(comp_fix_cost, f) << "," << diff(pay_fix_cost, f) << "," 
		<< diff(comp_prem_rec, f) << "," << diff(pay_prem_rec, f) << "," << diff(comp_prof, f) << "," << diff(pay_prof, f) << "," << diff(comp_act_prof, f) << "," << diff(pay_act_prof, f) << "," << diff(avg_brd_util, f) << "," << diff(avg_util, f) << endl;
	}
	file.close();

	string out_file = out_path + "results_full.csv";
	file.open(out_file);
	file << "result,num_prf,num_pat,num_prv,num_plnc,num_pln,num_self_insur,num_comp_insur,num_pay_insur,num_tot_insur,num_pln_off_comp,num_pln_off,pat_tot_pay,pat_tot_cost,self_pat_cost,";
	file << "comp_pat_cost,pay_pat_cost,prv_tot_cont,comp_prv_cont,pay_prv_cont,tot_prv_fix_cost,comp_fix_cost,pay_fix_cost,comp_prem_rec,pay_prem_rec,comp_prof,pay_prof,act_comp_prof,act_pay_prof,avg_brd_util,avg_util" << endl; 
	for (int f = 0; f < 3; ++f) {
		if (f == 0)
			file << "single,";
		else if (f == 1)
			file << "pre-opt,";
		else
			file << "post-opt,";

		file << fixed << setprecision(0) << num_pat << "," << tot_num_pat << "," << num_prv << "," << num_plnc << "," << num_pln << ","
		<< num_self_insur[f] << "," << num_comp_insur[f] << "," << num_pay_insur[f] << "," << num_tot_insur[f] << "," << num_pln_off_comp[f] << "," << num_pln_off[f] << "," << pat_tot_pay[f] << "," << pat_tot_cost[f] 
		<< "," << self_pat_cost[f] << "," << comp_pat_cost[f] << "," << pay_pat_cost[f] << "," << prv_tot_cont[f] << "," << comp_prv_cont[f] << "," << pay_prv_cont[f] << "," 
		<< prv_fix_cost[f] << "," << comp_fix_cost[f] << "," << pay_fix_cost[f] << "," << comp_prem_rec[f] << "," << pay_prem_rec[f] << "," << comp_prof[f] << "," << pay_prof[f] << ","
		<< comp_act_prof[f] << "," << pay_act_prof[f] << "," << avg_brd_util[f] << "," << avg_util[f] << endl;
	}
	file.close();

}

double diff(vector<int>& vec, int pos1, int denom) {
	
	if (denom == 0)
		denom = vec[pos1];
	if (denom == 0)
		denom = 1;

	return round_p((vec[pos1 + 1] - vec[pos1]) * 100.0 / abs(denom), 2);
}

double diff(vector<int>& vec, int pos1) {
	
	return vec[pos1 + 1] - vec[pos1];
}

double diff(vector<double>& vec, int pos1, int denom) {
	
	if (denom == 0)
		denom = vec[pos1];
	if (denom == 0)
		denom = 1;

	return round_p((vec[pos1 + 1] - vec[pos1]) * 100.0 / abs(denom), 2);
}

double diff(vector<double>& vec, int pos1) {
	
	return vec[pos1 + 1] - vec[pos1];
}








