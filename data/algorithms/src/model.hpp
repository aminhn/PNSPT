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

#include "header.hpp"
#include <ilcplex/ilocplex.h>

void build_MP(IloCplex& MP_cpl);

struct MP_var_struct {

	IloNumVarArray xj;
	vector<IloNumVarArray> zip;

	MP_var_struct(IloEnv env) {
		
		xj = IloNumVarArray(env, num_prv, 0, 1, ILOBOOL);

		zip.reserve(num_pat);
		for (int i = 0; i < num_pat; ++i) 
			zip.emplace_back(env, num_pln - 1, 0, 1, ILOBOOL);
	}

	MP_var_struct() {};

};


extern MP_var_struct MP_vars;
