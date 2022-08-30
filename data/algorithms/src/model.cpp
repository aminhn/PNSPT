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

#include "model.hpp"

MP_var_struct MP_vars;

void build_MP(IloCplex& MP_cpl) {

	IloEnv env = MP_cpl.getEnv();
	IloModel MP_mod = MP_cpl.getModel();

	MP_vars = MP_var_struct(env);

	IloExpr obj_expr(env);

	for (int i = 0; i < num_prf; ++i) {
		for (int p = 1; p < num_pln; ++p) 
			obj_expr += (plans[p]->prem[i] - plans[p]->cost[i]) * profiles[i]->pop * MP_vars.zip[i][p - 1];
	}

	for (int j = 0; j < num_prv; ++j)
		obj_expr -= providers[j]->ann_cost * MP_vars.xj[j];

	MP_mod.add(IloMaximize(env, obj_expr));
	obj_expr.end();


	for (int i = 0; i < num_prf; ++i) {
		for (int p = 1; p < num_pln; ++p) {
			for (vector<int>::iterator it = plans[p]->providers.begin(); it != plans[p]->providers.end(); ++it) 
				MP_mod.add(MP_vars.zip[i][p - 1] - MP_vars.xj[*it] <= 0);
		}
	}

	for (int i = 0; i < num_pat; ++i) {
		IloExpr temp_expr(env);
		for (int p = 1; p < num_pln; ++p) 
			temp_expr += MP_vars.zip[i][p - 1];
		MP_mod.add(temp_expr <= 1);
		temp_expr.end();
	}


	/*for (int i = 0; i < num_prf; ++i) {
		for (int p = 1; p < num_pln; ++p) {
			if (plans[p]->prem_comp[i] == profiles[i]->price_out) 
				MP_mod.add(MP_vars.zip[i][p - 1] <= 0);
		}
	}*/


}
