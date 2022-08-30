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
#include <iomanip>

void rec_data(IloCplex& model, bool fin);

double pre_obj = 0, impr_obj = 0, min_prec = 50000, org_obj;
int sol_count = 0;

void col_gen() {

	IloEnv MP_env;
	IloModel MP_mod(MP_env);
	IloCplex MP_cpl(MP_mod);

	IloEnv SP_env;
	IloModel SP_mod(SP_env);
	IloCplex SP_cpl(SP_mod);

	build_MP(MP_cpl);
	MP_cpl.setOut(MP_cpl.getEnv().getNullStream());

	cout << fixed << setprecision(0);

	rec_data(MP_cpl, 0);

	for(int iter = 0; ; ++iter) {

		cout << "Solving Master Problem\n";
		MP_cpl.solve();
		cout << "Done\n\n";

		if (MP_cpl.getStatus() != IloAlgorithm::Optimal) {
			cout << "MP " << MP_cpl.getStatus();
			return;
		}

		if (iter > 0 && floor(impr_obj + pre_obj) != floor(MP_cpl.getObjValue())) {
			cout << "\n\n!!!!!!!!!!!!!!!!!!!!!! Improvement in Objective Does Not Match Subproblem !!!!!!!!!!!!!!!!!!!!\n\n";
			cout << "Prev Obj: " << pre_obj << " " << "  SubP Obj: " << impr_obj << "  Sum: " << impr_obj + pre_obj << "  New Obj: " << MP_cpl.getObjValue() << endl;
		}

		pre_obj = MP_cpl.getObjValue();

		MP_vars.ass_vals(MP_cpl);


		ass_shdprc(MP_cpl);

		if (iter == 0) {
			org_obj = pre_obj;
			
			if (MP_cpl.getObjValue() > 0)
				min_prec = round(0.001 * MP_cpl.getObjValue());

			IloModel SP2_mod(SP_env);
			IloCplex SP2_cpl(SP2_mod);
			SP2_cpl.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, min_prec);
			SP2_cpl.setParam(IloCplex::Param::TimeLimit, 3600);
			SP2_cpl.setParam(SP2_cpl.Threads, 4);

			cout << "building Sub Problem 2\n";
			build_SP2(SP2_cpl, MP_cpl);
						
			cout << "Solving Sub Problem 2\n";
			SP2_cpl.solve();

			cout << "Done\n\n";
			if (SP2_cpl.getStatus() == IloAlgorithm::Optimal || SP2_cpl.getStatus() == IloAlgorithm::Feasible) {	
				SP_vars.ass_vals(SP2_cpl, 0);
				impr_obj = SP2_cpl.getObjValue();
			}
			else {
				cout << "SP2 " << SP2_cpl.getStatus() << endl;
				impr_obj = -1;
			}
		}
		else if (iter == 1) {
			
			SP_cpl.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, min_prec);
			SP_cpl.setParam(SP_cpl.Threads, 4);
			SP_cpl.setParam(IloCplex::Param::TimeLimit, 3600);
			
			cout << "building Sub Problem\n";
			build_SP(SP_cpl);
			
			cout << "Solving Sub Problem\n";
			SP_cpl.solve();
			cout << "Done\n\n";
			
			if (SP_cpl.getStatus() == IloAlgorithm::Optimal || SP_cpl.getStatus() == IloAlgorithm::Feasible) {
				SP_vars.ass_vals(SP_cpl, sol_count);
				++sol_count;
				impr_obj = SP_cpl.getObjValue();
			}
			else {
				cout << "SP " << SP_cpl.getStatus() << endl;
				impr_obj = 0;
			}
		}
		else {
			if (sol_count == SP_cpl.getNMIPStarts())
				impr_obj = 0;

			while (sol_count < SP_cpl.getNMIPStarts()) {
				SP_vars.ass_vals(SP_cpl, sol_count);
				++sol_count;

				impr_obj = cal_SP_obj();

				if (impr_obj >= min_prec) {
					cout << "Reduced cost column found from another SP solution!\n";
					break;
				}
			}
		
			if (sol_count >= SP_cpl.getNMIPStarts() && impr_obj < min_prec) {
				
				modify_SP(SP_cpl);
				
				cout << "Solving Sub Problem\n";
				SP_cpl.solve();
				cout << "Done\n\n";
				
				if (SP_cpl.getStatus() == IloAlgorithm::Optimal || SP_cpl.getStatus() == IloAlgorithm::Feasible) {
					SP_vars.ass_vals(SP_cpl, 0);
					sol_count = 1;
					impr_obj = SP_cpl.getObjValue();
				}
				else {
					cout << "SP " << SP_cpl.getStatus() << endl;
					impr_obj = 0;
				}
			}
		}

		if (impr_obj >= min_prec) {
			cout << "Positive reduced cost column exists with objective: " << impr_obj << endl;
			cout << "Improvement in objective over single plan: " << abs((pre_obj + impr_obj - max_single_prof) / max_single_prof * 100) << "%" << endl;
			cout << "Actual Improvement in objective over single plan: " << abs((pre_obj + impr_obj - max_single_prof) / (max_single_prof - deduct) * 100) << "%" << endl;
			cout << "Improvement in objective: " << abs((pre_obj + impr_obj - org_obj) / org_obj * 100) << "%" << endl;
			cout << "Actual Improvement in objective: " << abs((pre_obj + impr_obj - org_obj) / (org_obj - deduct) * 100) << "%" << endl;
			if (!update_plans()) {
				cout << "***************************Column already exists in MP!****************************\n";
				
				//rec_data(MP_cpl, 1);
				//return;
			}

			add_column(MP_cpl);
			SP_vars.hp.add(IloNumVar(SP_env, 0, 1, ILOBOOL));

		}
		else {
			if (impr_obj != -1) {
				cout << "objective: " << impr_obj << endl;
				cout << "**************************No Positive reduced cost column exists from SP1!!*********************************\n";

				cout << "Building SP2...\n";
				IloModel SP2_mod(SP_env);
				IloCplex SP2_cpl(SP2_mod);
				SP2_cpl.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, min_prec);
				SP2_cpl.setParam(IloCplex::Param::TimeLimit, 3600);
				SP2_cpl.setParam(SP2_cpl.Threads, 4);

				build_SP2(SP2_cpl, MP_cpl);

				cout << "Solving Sub Problem 2\n";
				SP2_cpl.solve();

				cout << "Done\n\n";
				if (SP2_cpl.getStatus() == IloAlgorithm::Optimal || SP2_cpl.getStatus() == IloAlgorithm::Feasible) {	
					SP_vars.ass_vals(SP2_cpl, 0);
					impr_obj = SP2_cpl.getObjValue();
				}
				else
					cout << "SP2 " << SP2_cpl.getStatus() << endl;
			}

			if (impr_obj >= min_prec) {
				cout << "Positive reduced cost column exists with objective: " << impr_obj << endl;
				cout << "Improvement in objective over single plan: " << abs((pre_obj + impr_obj - max_single_prof) / max_single_prof * 100) << "%" << endl;
				cout << "Actual Improvement in objective over single plan: " << abs((pre_obj + impr_obj - max_single_prof) / (max_single_prof - deduct) * 100) << "%" << endl;
				cout << "Improvement in objective: " << abs((pre_obj + impr_obj - org_obj) / org_obj * 100) << "%" << endl;
				cout << "Act Improvement in objective: " << abs((pre_obj + impr_obj - org_obj) / (org_obj - deduct) * 100) << "%" << endl;
				if (!update_plans()) 
					cout << "***************************Column already exists in MP!****************************\n";

				sol_count = SP_cpl.getNMIPStarts();
				add_column(MP_cpl);
				SP_vars.hp.add(IloNumVar(SP_env, 0, 1, ILOBOOL));
			}
			else {
				cout << "**************************No Positive reduced cost column exists from SP2!!*********************************\n";
				cout << "Improvement in objective over single plan: " << abs((pre_obj + impr_obj - max_single_prof) / max_single_prof * 100) << "%" << endl;
				cout << "Actual Improvement in objective over single plan: " << abs((pre_obj + impr_obj - max_single_prof) / (max_single_prof - deduct) * 100) << "%" << endl;
				cout << "Improvement in objective: " << ((pre_obj + impr_obj - org_obj) / org_obj * 100) << "%" << endl;
				cout << "Act Improvement in objective: " << abs((pre_obj + impr_obj - org_obj) / (org_obj - deduct) * 100) << "%" << endl;
				rec_data(MP_cpl, 1);
				return;
			}
		}
	
	}

}

void rec_data(IloCplex& MP_cpl, bool fin) {

	IloModel IP_model(MP_cpl.getEnv());
	IP_model.add(MP_cpl.getModel());

	for (int i = 0; i < num_pat; ++i) {
		for (int p = 0; p < num_plnc; ++p) 
			IP_model.add(IloConversion(IP_model.getEnv(), MP_vars.yip_comp[i][p], ILOBOOL));
		for (int p = 0; p < num_pln; ++p) 
			IP_model.add(IloConversion(IP_model.getEnv(), MP_vars.yip[i][p], ILOBOOL));
	}

	IloCplex IP_cpl(IP_model);
	IP_cpl.setOut(IP_cpl.getEnv().getNullStream());
	IP_cpl.solve();

	cout << "**********IP MODEL***********\n";
	if (IP_cpl.getStatus() != IloAlgorithm::Optimal) {
		cout << "IP " << IP_cpl.getStatus();
		return;
	}

	cout << "Objective: " << IP_cpl.getObjValue() << endl;

	rec_results(IP_cpl, fin);
	IP_cpl.end();
	IP_model.end();
}


