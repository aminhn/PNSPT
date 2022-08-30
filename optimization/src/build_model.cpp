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
#include <algorithm>

vector<IloRange> sum_con;
vector<vector<vector<IloRange>*>*> linkzx;
vector<vector<IloRange>*> util_conp;
vector<vector<IloRange>*> util_conpc;

vector<IloRange> zip_sum;
vector<vector<IloRange>*> linkyz;

vector<vector<bool>> utilc_indic;
bool constructed = 0;

vector<bool> target_pat;

MP_var_struct MP_vars;
SP_var_struct SP_vars;

vector<double> fjxj;
vector<double> riciyi;

void build_MP(IloCplex& MP_cpl) {

	IloEnv env = MP_cpl.getEnv();
	IloModel MP_mod = MP_cpl.getModel();

	MP_vars = MP_var_struct(env);

	IloExpr obj_expr(env);

	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_pln; ++p) 
			obj_expr += (plans[p]->prem[i] - plans[p]->cost[i]) * patients[i]->pop * MP_vars.yip[i][p];
	}

	for (int j = 0; j < num_prv; ++j)
		obj_expr -= providers[j]->ann_cost * MP_vars.xj[j];

	MP_mod.add(IloMaximize(env, obj_expr));
	obj_expr.end();


	zip_sum.reserve(num_pat);
	for (int i = 0; i < num_pat; ++i) {
		IloExpr temp_expr(env);
		for (int p = 0; p < num_pln; ++p) 
			temp_expr += MP_vars.zip[i][p];
		zip_sum.emplace_back(temp_expr <= 1);
		MP_mod.add(zip_sum.back());
		temp_expr.end();
	}


	linkzx.reserve(num_pat);
	for (int i = 0; i < num_pat; ++i) {
		linkzx.emplace_back(new vector<vector<IloRange>*>);
		linkzx.back()->reserve(num_pln);
		for (int p = 0; p < num_pln; ++p) {
			linkzx.back()->emplace_back(new vector<IloRange>);
			linkzx.back()->back()->reserve(plans[p]->providers.size());
			for (vector<int>::iterator it = plans[p]->providers.begin(); it != plans[p]->providers.end(); ++it) {
				linkzx.back()->back()->emplace_back(MP_vars.zip[i][p] - MP_vars.xj[*it] <= 0);
				MP_mod.add(linkzx.back()->back()->back());
			}
		}
	}

	linkyz.reserve(num_pat);
	for (int i = 0; i < num_pat; ++i) {
		linkyz.emplace_back(new vector<IloRange>);
		linkyz.back()->reserve(num_pln);
		for (int p = 0; p < num_pln; ++p) {
			linkyz[i]->emplace_back(MP_vars.yip[i][p] - MP_vars.zip[i][p] <= 0);
			MP_mod.add(linkyz[i]->back());
		}
	}


	sum_con.reserve(num_pat);
	for (int i = 0; i < num_pat; ++i) {
		IloExpr temp_expr(env);
		for (int p = 0; p < num_plnc; ++p)
			temp_expr += MP_vars.yip_comp[i][p];
		for (int p = 0; p < num_pln; ++p)
			temp_expr += MP_vars.yip[i][p];

		sum_con.emplace_back(temp_expr == 1);
		MP_mod.add(sum_con.back());
		temp_expr.end();
	}


	util_conp.reserve(num_pat);
	for (int i = 0; i < num_pat; ++i) {
		util_conp.emplace_back(new vector<IloRange>);
		util_conp[i]->reserve(num_pln);
		for (int p = 0; p < num_pln; ++p) {
			IloExpr temp_expr(env);
			temp_expr += MP_vars.zip[i][p];
			for (int p2 = 0; p2 < num_plnc; ++p2) {
				if (comp_plans[p2]->util[i] + 1 < plans[p]->util[i]) 
					temp_expr += MP_vars.yip_comp[i][p2];
			}
			for (int p2 = 0; p2 < num_pln; ++p2) {
				if (plans[p2]->util[i] + 1 < plans[p]->util[i]) 
					temp_expr += MP_vars.yip[i][p2];
			}
			util_conp[i]->emplace_back(temp_expr <= 1);
			MP_mod.add(util_conp[i]->back());
			temp_expr.end();
		}
	}

	util_conpc.reserve(num_pat);
	utilc_indic = vector<vector<bool>>(num_pat, vector<bool>(num_plnc, 1));
	for (int i = 0; i < num_pat; ++i) {
		util_conpc.emplace_back(new vector<IloRange>);
		util_conpc[i]->reserve(num_plnc);
		for (int p = 0; p < num_plnc; ++p) {
			IloExpr temp_expr(env);
			bool empty = 1;
			for (int p2 = 0; p2 < num_plnc; ++p2) {
				if (comp_plans[p2]->util[i] + 1 < comp_plans[p]->util[i]) {
					temp_expr += MP_vars.yip_comp[i][p2];
					empty = 0;
				}
			}
			for (int p2 = 0; p2 < num_pln; ++p2) {
				if (plans[p2]->util[i] + 1 < comp_plans[p]->util[i]) {
					temp_expr += MP_vars.yip[i][p2];
					empty = 0;
				}
			}

			if (!empty) {
				util_conpc[i]->emplace_back(temp_expr <= 0);
				MP_mod.add(util_conpc[i]->back());
			}
			else {
				util_conpc[i]->emplace_back(IloRange());
				utilc_indic[i][p] = 0;
			}
			temp_expr.end();
		}
	}

}


void build_SP(IloCplex& SP_cpl) {

	IloEnv env = SP_cpl.getEnv();
	IloModel SP_mod = SP_cpl.getModel();

	if (!constructed)
		SP_vars = SP_var_struct(env);

	IloExpr obj_expr(env);
	for (int i = 0; i < num_pat; ++i) {
		obj_expr += patients[i]->pop * SP_vars.ri[i] - patients[i]->pop * riciyi[i] * SP_vars.zi[i];

		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser; ++s)
				obj_expr -= patients[i]->pop * patients[i]->cost_prv[j]->at(s) * SP_vars.wijs[i][j][s];
		}
	}

	for (int j = 0; j < num_prv; ++j) 
		obj_expr -= (providers[j]->ann_cost - fjxj[j]) * SP_vars.vj[j];

	SP_mod.add(IloMaximize(env, obj_expr));

	obj_expr.end();

	for (int i = 0; i < num_pat; ++i) {
		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser; ++s) {
				if (providers[j]->services[s])
					SP_mod.add(SP_vars.wijs[i][j][s] - SP_vars.vj[j] <= 0);
				else
					SP_mod.add(SP_vars.wijs[i][j][s] <= 0);
			}
		}
	}

	for (int i = 0; i < num_pat; ++i) {
		for (int s = 0; s < num_ser; ++s) {
			IloExpr temp_expr(env);
			temp_expr += SP_vars.wi0s[i][s] - SP_vars.zi[i];
			for (int j = 0; j < num_prv; ++j)
				temp_expr += SP_vars.wijs[i][j][s];

			SP_mod.add(temp_expr == 0);
			temp_expr.end();
		}
	}


	for (int i = 0; i < num_pat; ++i) {
		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser; ++s) {
				if (providers[j]->services[s]) {
					IloExpr temp_expr(env);
					temp_expr += SP_vars.wi0s[i][s] + SP_vars.vj[j];
					for (int jp = 0; jp < num_prv; ++jp) {
						if (providers[jp]->services[s] && (patients[i]->util_prv_ser[jp]->at(s) < patients[i]->util_prv_ser[j]->at(s) || (patients[i]->util_prv_ser[jp]->at(s) == patients[i]->util_prv_ser[j]->at(s) && patients[i]->cost_prv[jp]->at(s) > patients[i]->cost_prv[j]->at(s)))) 
							temp_expr += SP_vars.wijs[i][jp][s];
					}
					SP_mod.add(temp_expr <= 1);
					temp_expr.end();
				}
			}
		}
	}

	for (int i = 0; i < num_pat; ++i)
		SP_mod.add(SP_vars.ri[i] - min(12000, patients[i]->price_out) * SP_vars.zi[i] <= 0);


	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_plnc; ++p) {
			double disc = 100;
			if (p == 0)
				disc = epsilon;
			IloExpr temp_expr(env);
			temp_expr -= patients[i]->beta_prc * SP_vars.ri[i] + (comp_plans[p]->util[i] + patients[i]->beta_prc * disc) * SP_vars.zi[i];
			for (int s = 0; s < num_ser; ++s)
				temp_expr -= patients[i]->beta_prc * patients[i]->cost_prv_self[patients[i]->glob_ch[s]]->at(s) * SP_vars.wi0s[i][s];
			for (int j = 0; j < num_prv; ++j) 
				temp_expr += patients[i]->beta_brd * patients[i]->util_prv[j] * SP_vars.vj[j];
			SP_mod.add(temp_expr >= 0);
			temp_expr.end();
		}
	}

}


void build_SP2(IloCplex& SP2_cpl, IloCplex& MP_cpl) {

	IloEnv env = SP2_cpl.getEnv();
	IloModel SP2_mod = SP2_cpl.getModel();

	if (!constructed)
		SP_vars = SP_var_struct(env);

	IloExpr obj_expr(env);
	for (int i = 0; i < num_pat; ++i) {
		obj_expr += patients[i]->pop * SP_vars.ri[i] - patients[i]->pop * riciyi[i] * SP_vars.zi[i];

		for (int p = 0; p < num_pln; ++p) {
			if (MP_cpl.getValue(MP_vars.zip[i][p]) > 0.1) 
				obj_expr += patients[i]->pop * riciyi[i] * SP_vars.zi[i] * SP_vars.hp[p] - patients[i]->pop * riciyi[i] * SP_vars.hp[p];
		}

		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser; ++s)
				obj_expr -= patients[i]->pop * patients[i]->cost_prv[j]->at(s) * SP_vars.wijs[i][j][s];
		}
	}

	for (int j = 0; j < num_prv; ++j) 
		obj_expr += fjxj[j] * SP_vars.aj[j] - (providers[j]->ann_cost - fjxj[j]) * SP_vars.vj[j];

	SP2_mod.add(IloMaximize(env, obj_expr));
	obj_expr.end();

	for (int i = 0; i < num_pat; ++i) {
		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser; ++s) {
				if (providers[j]->services[s])
					SP2_mod.add(SP_vars.wijs[i][j][s] - SP_vars.vj[j] <= 0);
				else
					SP2_mod.add(SP_vars.wijs[i][j][s] <= 0);
			}
		}
	}

	for (int i = 0; i < num_pat; ++i) {
		for (int s = 0; s < num_ser; ++s) {
			IloExpr temp_expr(env);
			temp_expr += SP_vars.wi0s[i][s] - SP_vars.zi[i];
			for (int j = 0; j < num_prv; ++j)
				temp_expr += SP_vars.wijs[i][j][s];

			SP2_mod.add(temp_expr == 0);
			temp_expr.end();
		}
	}

	for (int i = 0; i < num_pat; ++i) {
		for (int j = 0; j < num_prv; ++j) {
			for (int s = 0; s < num_ser; ++s) {
				if (providers[j]->services[s]) {
					IloExpr temp_expr(env);
					temp_expr += SP_vars.wi0s[i][s] + SP_vars.vj[j];
					for (int jp = 0; jp < num_prv; ++jp) {
						if (providers[jp]->services[s] && (patients[i]->util_prv_ser[jp]->at(s) < patients[i]->util_prv_ser[j]->at(s) || (patients[i]->util_prv_ser[jp]->at(s) == patients[i]->util_prv_ser[j]->at(s) && patients[i]->cost_prv[jp]->at(s) > patients[i]->cost_prv[j]->at(s)))) 
							temp_expr += SP_vars.wijs[i][jp][s];
					}
					SP2_mod.add(temp_expr <= 1);
					temp_expr.end();
				}
			}
		}
	}

	for (int i = 0; i < num_pat; ++i)
		SP2_mod.add(SP_vars.ri[i] - min(12000, patients[i]->price_out) * SP_vars.zi[i] <= 0);

	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_plnc; ++p) {
			double disc = 100;
			if (p == 0)
				disc = epsilon;
			IloExpr temp_expr(env);
			temp_expr -= patients[i]->beta_prc * SP_vars.ri[i] + (comp_plans[p]->util[i] + patients[i]->beta_prc * disc) * SP_vars.zi[i];
			for (int s = 0; s < num_ser; ++s)
				temp_expr -= patients[i]->beta_prc * patients[i]->cost_prv_self[patients[i]->glob_ch[s]]->at(s) * SP_vars.wi0s[i][s];
			for (int j = 0; j < num_prv; ++j)
				temp_expr += patients[i]->beta_brd * patients[i]->util_prv[j] * SP_vars.vj[j];
			SP2_mod.add(temp_expr >= 0);
			temp_expr.end();
		}
	}

	for (int j = 0; j < num_prv; ++j) 
		SP2_mod.add(SP_vars.aj[j] + SP_vars.vj[j] <= 1);

	for (int p = 0; p < num_pln; ++p) {
		for (vector<int>::iterator it = plans[p]->providers.begin(); it != plans[p]->providers.end(); ++it) 
			SP2_mod.add(SP_vars.aj[*it] - SP_vars.hp[p] <= 0);
	}
}

void ass_shdprc(IloCplex& MP_cpl) {

	cout << MP_cpl.getStatus() << ", Obj value: " << MP_cpl.getObjValue() << endl;

	cout << " xj:\n";
	for (int j = 0; j < num_prv; ++j) 
		cout << abs(MP_vars.xj_val[j]) << " ";
	
	cout << endl;
	cout << "     yip_comp:";
	for (int p = 0; p < num_plnc - 1; ++p)
		cout << "  ";
	cout << "zip:   ";
	for (int p = 0; p < num_pln; ++p)
		cout << "  ";
	/*cout << "yip:   ";
	for (int p = 0; p < num_pln; ++p)
		cout << "  ";
	cout << "uip:";*/
	cout << endl;
	for (int i = 0; i < num_pat; ++i) {
		cout << i + 1 << " ";
		for (int pp = 0; pp <= floor(log10(num_pat) + 1) - floor(log10(i + 1) + 1); ++pp)
			cout << " ";
		for (int p = 0; p < num_plnc; ++p) 
			cout << abs(MP_vars.yip_comp_val[i][p]) << " ";
		cout << "       ";
		for (int p = 0; p < num_pln; ++p) 
			cout << abs(MP_vars.zip_val[i][p]) << " ";
		cout << "       ";
		for (int p = 0; p < num_pln; ++p) 
			cout << abs(MP_vars.yip_val[i][p]) << " ";
		cout << "       ";
		/*for (int p = 0; p < num_plnc; ++p) 
			cout << comp_plans[p]->util[i] << " ";
		cout << "       ";
		for (int p = 0; p < num_pln; ++p) 
			cout << plans[p]->util[i] << " ";
		*/
		cout << endl;
	}

	double tot_income = 0, tot_cost = 0, tot_fix_cost = 0;
	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_pln; ++p) {
			if (MP_vars.yip_val[i][p] > 0.5) {
				tot_income += plans[p]->prem[i] * patients[i]->pop;
				tot_cost += plans[p]->cost[i] * patients[i]->pop;
			}
		}
	}
	for (int j = 0; j < num_prv; ++j) {
		if (MP_vars.xj_val[j] > 0.5) 
			tot_fix_cost += providers[j]->ann_cost;
	}

	cout << "Total income: " << tot_income << "\nTotal cost:" << tot_cost << "\nTotal fixed cost: " << tot_fix_cost << endl;

	for (int i = 0; i < num_pat; ++i) {
		if (!target_pat.empty() && MP_vars.yip_val[i][num_pln - 1] != target_pat[i]) {
			cout << "\n!!!!! Targeted patient " << i + 1 << " picked wrong plan!!!!!\n";
			cin.get();
		}
	}


	fjxj = vector<double>(num_prv, 0);
	
	for (int j = 0; j < num_prv; ++j) {
		if (MP_vars.xj_val[j] > 0.1)
			fjxj[j] = providers[j]->ann_cost;
	}

	riciyi = vector<double>(num_pat, 0);

	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_pln; ++p) {
			if (MP_vars.yip_val[i][p] > 0.1) {
				riciyi[i] = plans[p]->prem[i] - plans[p]->cost[i];
				break;
			}
		}
	}

}

void add_column(IloCplex& MP_cpl) {

	IloEnv MP_env = MP_cpl.getEnv();
	IloModel MP_mod = MP_cpl.getModel();
	IloObjective MP_obj = MP_cpl.getObjective();

	cout << "New Plan: \n";
	for (int j = 0; j < num_prv; ++j) {
		cout << abs(SP_vars.vj_val[j]) << " ";
		//if (SP_vars.vj_val[j] > 0.5)
			//cout << j << " ";
	}
	cout << endl;

	target_pat = vector<bool>(num_pat, 0);
	cout << "Targetted patients: \n";
	for (int i = 0; i < num_pat; ++i) {
		if (SP_vars.zi_val[i] > 0.5 && plans.back()->prem[i] - plans.back()->cost[i] > 0) {
			target_pat[i] = 1;
			cout << i << " ";
			for (int p = 0; p < num_pln - 1; ++p)
				MP_obj.setLinearCoef(MP_vars.yip[i][p], 0);
		}
	}
	cout << endl;

	for (int i = 0; i < num_pat; ++i) {
		MP_vars.zip[i].add(IloNumVar(MP_env, 0, 1, ILOFLOAT));
		MP_vars.yip[i].add(IloNumVar(MP_env, 0, 1, ILOFLOAT));
	}

	for (int i = 0; i < num_pat; ++i) 
		MP_obj.setLinearCoef(MP_vars.yip[i][num_pln - 1], (plans.back()->prem[i] - plans.back()->cost[i]) * patients[i]->pop);

	for (int i = 0; i < num_pat; ++i)
		zip_sum[i].setExpr(zip_sum[i].getExpr() + MP_vars.zip[i][num_pln - 1]);
	
	for (int i = 0; i < num_pat; ++i) {
		linkzx[i]->emplace_back(new vector<IloRange>);
		linkzx[i]->back()->reserve(plans.back()->providers.size());
		for (vector<int>::iterator it = plans.back()->providers.begin(); it != plans.back()->providers.end(); ++it) {
			linkzx[i]->back()->emplace_back(MP_vars.zip[i][num_pln - 1] - MP_vars.xj[*it] <= 0);
			MP_mod.add(linkzx[i]->back()->back());
		}
	}

	for (int i = 0; i < num_pat; ++i) {
		linkyz[i]->emplace_back(MP_vars.yip[i][num_pln - 1] - MP_vars.zip[i][num_pln - 1] <= 0);
		MP_mod.add(linkyz[i]->back());
	}
	
	for (int i = 0; i < num_pat; ++i)
		sum_con[i].setExpr(sum_con[i].getExpr() + MP_vars.yip[i][num_pln - 1]);

	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_pln - 1; ++p) {
			if (plans.back()->util[i] + 1 < plans[p]->util[i]) 
				util_conp[i]->at(p).setExpr(util_conp[i]->at(p).getExpr() + MP_vars.yip[i][num_pln - 1]);
		}
		IloExpr temp_expr(MP_env);
		temp_expr += MP_vars.zip[i][num_pln - 1];
		for (int p = 0; p < num_plnc; ++p) {
			if (comp_plans[p]->util[i] + 1 < plans.back()->util[i]) 
				temp_expr += MP_vars.yip_comp[i][p];
		}
		for (int p = 0; p < num_pln - 1; ++p) {
			if (plans[p]->util[i] + 1 < plans.back()->util[i]) 
				temp_expr += MP_vars.yip[i][p];
		}
		util_conp[i]->emplace_back(temp_expr <= 1);
		MP_mod.add(util_conp[i]->back());
		temp_expr.end();
	}

	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_plnc; ++p) {
			if (plans.back()->util[i] + 1 < comp_plans[p]->util[i]) {
				if (utilc_indic[i][p])
					util_conpc[i]->at(p).setExpr(util_conpc[i]->at(p).getExpr() + MP_vars.yip[i][num_pln - 1]);
				else {
					util_conpc[i]->at(p) = MP_vars.yip[i][num_pln - 1] <= 0;
					MP_mod.add(util_conpc[i]->at(p));
					utilc_indic[i][p] = 1;
				}
			}
		}
	}

}

void modify_SP(IloCplex& SP_cpl) {

	IloObjective SP_obj = SP_cpl.getObjective();

	for (int j = 0; j < num_prv; ++j) 
		SP_obj.setLinearCoef(SP_vars.vj[j], fjxj[j] - providers[j]->ann_cost);

	for (int i = 0; i < num_pat; ++i) 
		SP_obj.setLinearCoef(SP_vars.zi[i], -patients[i]->pop * riciyi[i]);	


}

double cal_SP_obj() {

	double obj = 0;

	for (int i = 0; i < num_pat; ++i) {
		if (SP_vars.zi_val[i] > 0.5) {
			double prof = patients[i]->pop * (SP_vars.ri_val[i] - riciyi[i] * SP_vars.zi_val[i]);

			for (int j = 0; j < num_prv; ++j) {
				for (int s = 0; s < num_ser; ++s)
					prof -= patients[i]->pop * patients[i]->cost_prv[j]->at(s) * SP_vars.wijs_val[i][j][s];
			}

			if (prof > 0)
				obj += prof;
			else 
				SP_vars.zi_val[i] = 0;
		}
	}

	for (int j = 0; j < num_prv; ++j) 
		obj -= (providers[j]->ann_cost - fjxj[j]) * SP_vars.vj_val[j];

	return obj;

}

