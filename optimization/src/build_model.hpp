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

#include "load_inst.hpp"
#include <ilcplex/ilocplex.h>

void col_gen();
void build_MP(IloCplex& MP_cpl);
void build_SP(IloCplex& SP_cpl);
void build_SP2(IloCplex& SP2_cpl, IloCplex& MP_cpl);
void ass_shdprc(IloCplex& P_cpl);
void add_column(IloCplex& MP_cpl);
void modify_SP(IloCplex& SP_cpl);
double cal_SP_obj();


struct MP_var_struct {

	IloNumVarArray xj;
	vector<IloNumVarArray> zip;
	vector<IloNumVarArray> yip;
	vector<IloNumVarArray> yip_comp;

	IloNumArray xj_val;
	vector<IloNumArray> zip_val;
	vector<IloNumArray> yip_val;
	vector<IloNumArray> yip_comp_val;

	MP_var_struct(IloEnv env) {
		
		xj = IloNumVarArray(env, num_prv, 0, 1, ILOFLOAT);

		xj_val = IloNumArray(env);

		yip.reserve(num_pat);
		for (int i = 0; i < num_pat; ++i) {
			zip.emplace_back(env, num_pln, 0, 1, ILOFLOAT);
			yip.emplace_back(env, num_pln, 0, 1, ILOFLOAT);
			yip_comp.emplace_back(env, num_plnc, 0, 1, ILOFLOAT);
			yip_val.emplace_back(env);
			zip_val.emplace_back(env);
			yip_comp_val.emplace_back(env);
		}
	}

	MP_var_struct() {};

	void ass_vals(IloCplex& MP_cpl) {
		MP_cpl.getValues(xj, xj_val);
		for (int i = 0; i < num_pat; ++i) {
			MP_cpl.getValues(yip[i], yip_val[i]);
			MP_cpl.getValues(zip[i], zip_val[i]);
			MP_cpl.getValues(yip_comp[i], yip_comp_val[i]);
		}
	}
};

struct SP_var_struct {

	IloNumVarArray ri;
	IloNumVarArray vj;
	IloNumVarArray zi;
	IloNumVarArray aj;
	IloNumVarArray hp;
	vector<IloNumVarArray> wi0s;
	vector<vector<IloNumVarArray>> wijs;

	IloNumArray ri_val;
	IloNumArray vj_val;
	IloNumArray zi_val;
	IloNumArray aj_val;
	IloNumArray hp_val;
	vector<IloNumArray> wi0s_val;
	vector<vector<IloNumArray>> wijs_val;

	SP_var_struct(IloEnv env) {

		ri = IloNumVarArray(env, num_pat, 0, IloInfinity, ILOFLOAT);
		vj = IloNumVarArray(env, num_prv, 0, 1, ILOBOOL);
		zi = IloNumVarArray(env, num_pat, 0, 1, ILOBOOL);
		aj = IloNumVarArray(env, num_prv, 0, 1, ILOBOOL);
		hp = IloNumVarArray(env, num_pln, 0, 1, ILOBOOL);

		ri_val = IloNumArray(env);
		vj_val = IloNumArray(env);
		zi_val = IloNumArray(env);
		aj_val = IloNumArray(env);
		hp_val = IloNumArray(env);


		wi0s.reserve(num_pat);
		for (int i = 0; i < num_pat; ++i) {
			wi0s.emplace_back(env, num_ser, 0, 1, ILOBOOL);
			wi0s_val.emplace_back(env);
		}


		wijs.reserve(num_pat);
		for (int i = 0; i < num_pat; ++i) {
			vector<IloNumVarArray> temp;
			vector<IloNumArray> temp2;
			wijs.push_back(temp);
			wijs_val.push_back(temp2);
			for (int j = 0; j < num_prv; ++j) {
				wijs[i].emplace_back(env, num_ser, 0, 1, ILOBOOL);
				wijs_val[i].emplace_back(env);
			}
		}
	}

	SP_var_struct() {};

	void ass_vals(IloCplex& SP, int sol_num) {
		SP.getValues(ri, ri_val, sol_num);
		SP.getValues(vj, vj_val, sol_num);
		SP.getValues(zi, zi_val, sol_num);
		for (int i = 0; i < num_pat; ++i)
			SP.getValues(wi0s[i], wi0s_val[i], sol_num);
		for (int i = 0; i < num_pat; ++i) {
			for (int j = 0; j < num_prv; ++j)
				SP.getValues(wijs[i][j], wijs_val[i][j], sol_num);
		}
	}

	void ass_vals2(IloCplex& SP, int sol_num) {
		SP.getValues(ri, ri_val, sol_num);
		SP.getValues(vj, vj_val, sol_num);
		SP.getValues(zi, zi_val, sol_num);
		SP.getValues(aj, aj_val, sol_num);
		SP.getValues(hp, hp_val, sol_num);
		for (int i = 0; i < num_pat; ++i)
			SP.getValues(wi0s[i], wi0s_val[i], sol_num);
		for (int i = 0; i < num_pat; ++i) {
			for (int j = 0; j < num_prv; ++j)
				SP.getValues(wijs[i][j], wijs_val[i][j], sol_num);
		}
	}


};


extern MP_var_struct MP_vars;
extern SP_var_struct SP_vars;




