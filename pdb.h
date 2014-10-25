#ifndef PDB_H
#define PDB_H

#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <map>

#include "string_operate.h"


struct  atom
{
	int atom_serial;
	char atom_name[20];
	char residue_name[10];
	int residue_serial;
	float x;
	float y;
	float z;
};



map<int, atom> read_pdb_to_atom(char * pdb_file_name)
{
	map<int, atom> atom_list;
	int n=0;
	string s;
	ifstream infile(pdb_file_name);
	while(!infile.fail())
	{
		getline(infile,s);
		if(s.find("ATOM")<100||s.find("HETATM")<100)
		{
			struct atom item;

			item.atom_serial=atoi(Split(s.substr(6,5),' ',0).c_str());
			strcpy(item.atom_name,Split(s.substr(12,4),' ',0).c_str());
			strcpy(item.residue_name,Split(s.substr(17,3),' ',0).c_str());
			item.residue_serial=atoi(Split(s.substr(22,4),' ',0).c_str());

			atom_list[item.atom_serial]= item;
		}
	}
	infile.close();
	return atom_list;
}

#endif


