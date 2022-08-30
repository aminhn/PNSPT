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

#include <time.h>
#include <string.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <time.h>
#include "load_inst.hpp"
#include "utility.hpp"
#include "build_model.hpp"

using namespace std;

string instance, in_path;

int main(int argc, char* argv[]) {

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] != '-')
			continue;
		else if (strcmp(argv[i], "-inst") == 0)
			instance = argv[i + 1];
	}

	in_path = instance + "/"; 

	load_inst();

	clock_t kk = clock();

	col_gen();

	cout << "\n\nAlgorithm Terminated Successfully\n\n";

	cout << "CPU time is: ";
	cout << "Total time: " << give_time(clock() - kk) << endl;

	return 0;


}


