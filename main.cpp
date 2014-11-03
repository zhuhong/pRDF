/*
first version 1.0.0, finished at 2014.10.25.
version 1.1.0. update 10.26.
	- delete pdb.h, choosinng solvent and solute from input.
	- renew the args.
*/



#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <utility>
#include <vector>
#include <map>
#include <iomanip>
#include <stdlib.h>

#include "string_operate.h"
#include "read_ndx.h"
// #include "pdb.h"
// #include "my_math.h"

using namespace std;

extern "C"
{
#include "xdrfile_xtc.h"
}


void Print_usage();



pair<float,int> Min_dist(rvec *x,vector<int> solu_l,float * grid, float * center)
{
	// '''
	// atom_l: A full atom list.
	// solu_l: A solute atom index list.
	// solv_i: Index for one solvent molecule.
	// '''
	float min_dist = 0.0;
	int   min_index =0;

	for (int i=0; i< solu_l.size();i++)
	{

		int item = solu_l[i];
		float tmp = sqrt(pow(x[item-1][0]-grid[0],2) \
				+ pow(x[item-1][1]-grid[1],2) \
				+ pow(x[item-1][2]-grid[2],2));
		if (i==0)
		{
			min_dist = tmp;
			min_index = item;
			continue;
		}
		if (tmp < min_dist)
		{
			min_dist = tmp;
			min_index = item;
		}  
	}
	pair<float,int> data(min_dist,min_index);

	return data;
}



vector< pair<float,int> > Sort_dist(rvec *x,vector<int> solu_l,float * grid,float * center)
{
	// vector<int> min_dist;
	// min_index =list()
	vector<pair<float,int> > dist_index;
	bool INSERT = false;
	for (int i=0; i< solu_l.size();i++)
	{
		int item = solu_l[i];
		float tmp = sqrt(pow(x[item-1][0]-grid[0],2) \
				+ pow(x[item-1][1]-grid[1],2) \
				+ pow(x[item-1][2]-grid[2],2));
		pair<float,int> aa(tmp,item);
		// cout<< tmp <<endl;
		if(i==0)
		{
			dist_index.push_back(aa);
			INSERT = true;
			continue;
		}
		int size = dist_index.size();
		INSERT = false;
		for(int i=0;i<size;i++)
		{
			if (tmp < dist_index[i].first)
			{
				dist_index.insert(dist_index.begin()+i,aa);
				// cout<<dist_index[i].first<<endl;
				INSERT = true;
				break;
			}
		}
		if (INSERT ==false)
		{
			dist_index.push_back(aa);
		}
	}
	return dist_index;
}



// vector<int> Get_solvent_list(map<int, atom> atom_list)
// {
// 	vector<int> solvent_list ;
// 	map<int,atom>::iterator it = atom_list.begin(); 
// 	for (;it!=atom_list.end(); ++it)
// 	{
// 		int i = it->first;
// 		// cout << it-> second.residue_name<< it->second.atom_name;
// 		if(strcmp(it->second.residue_name,"WAT") ==0 && strcmp(it->second.atom_name,"O")==0)
// 		{
// 			solvent_list.push_back(i);
// 		}
// 		else if(strcmp(it->second.residue_name,"SOL")==0 && strcmp(it->second.atom_name,"OW")==0)
// 		{
// 			solvent_list.push_back(i);
// 		}
// 	}
// 	return solvent_list;
// }

vector<float> hist(vector<float> value_list,vector<int> number_list,float min_value,float max_value,int nbins)
{
	float  bin       = (max_value - min_value) / nbins;
	// cout << bin << endl;
	vector<float> grid_list;
	vector<float> numb_list;
	for (int i=0;i< nbins;i++)
	{
		grid_list.push_back(0);
		numb_list.push_back(0);
	}
	// cout << value_list.size()<< endl;
	for(int i=0;i< value_list.size();i++)
	{
		int temp = (value_list[i] - min_value) / bin ;
		// cout<< temp<< endl;
		if (value_list[i] < max_value)
		{
			grid_list[temp] += 1;
			numb_list[temp] += number_list[i];
		}
	}
	for(int i=0;i<nbins;i++)
	{
		if (grid_list[i] > 0)
		{
			numb_list[i] /= (1.0*grid_list[i]);
			// cout<< numb_list[i]<< endl;
		}
	}
	return numb_list;
}

void pRDF(char * traj_file,char * coor_file,char * index_file, char * data_file, int dmax)
{

	// map<int,atom> atom_list    =read_pdb_to_atom(coor_file);
	int solute_index;
	int solvent_index;
	vector <Index_class> index_list   =Read_index_to_Inclass(index_file);
	Print_Index(index_list);
	cout << "Choosing solute:" ;
	cin >> solute_index; 
	vector<int> solute_list  =index_list[solute_index].group_list ;
	cout << "Choosing solvent: " ;
	cin >> solvent_index;
	vector<int> solvent_list = index_list[solvent_index].group_list ;
	// vector<int> solvent_list =Get_solvent_list(atom_list);
	// cout<< solvent_list.size()<< endl;

	float GR[100] ; //=numpy.zeros((100),dtype=numpy.float64)
	for(int i=0;i<100;i++)
	{
		GR[i] = 0.0;
	}
	int FRAMES = 0;
	ofstream Data_file(data_file);


	int natoms,step;
	float time_temp;
	float p;
	// vector<double> result;
	matrix box;
	rvec *x;
	XDRFILE *xtc;
	xtc=xdrfile_open(traj_file,"r");
	int read_return=read_xtc_natoms(traj_file,&natoms);
	x=(rvec * )calloc(natoms,sizeof(x[0]));
	while(1)
	{
		read_return=read_xtc(xtc,natoms,&step,&time_temp,box,x,&p);
		// if(step%100000==0)
		// {

		// 	cout<<"Reading frame"<<"\t"<<step<<" time "<<time_temp<<endl;
		// }
		if(read_return!=0)
		{
			break;
		}

		FRAMES += 1;
		cout<< "reading frame: " << FRAMES << endl;

		// solute_atoms = dict()
		float coor_center[3];
		coor_center[0]=0.0;
		coor_center[1]=0.0;
		coor_center[2]=0.0;
		float X_min = x[0][0]; 
		float X_max = x[0][0]; 
		float Y_min = x[0][1]; 
		float Y_max = x[0][1]; 
		float Z_min = x[0][2]; 
		float Z_max = x[0][2];         
		float R_solute =0.0;

		for (int ii=0;ii<solute_list.size();ii++)
		{
			int i = solute_list[ii];
			coor_center[0] += x[i-1][0];
			coor_center[1] += x[i-1][1];
			coor_center[2] += x[i-1][2];
			if (x[i-1][0] < X_min)
			{
				X_min = x[i-1][0];
			}
			else if (x[i-1][0] > X_max)
			{
				X_max = x[i-1][0];
			}

			if (x[i-1][1] < Y_min)
			{
				Y_min = x[i-1][1];
			}
			else if (x[i-1][1] > Y_max)
			{
				Y_max = x[i-1][1];
			}

			if (x[i-1][2] < Z_min)
			{
				Z_min = x[i-1][2];
			}
			else if (x[i-1][2] > Z_max)
			{
				Z_max =x[i-1][2];
			}
		}

		coor_center[0] /=solute_list.size();
		coor_center[1] /=solute_list.size(); 
		coor_center[2] /=solute_list.size();


		for (int ii=0;ii<solute_list.size();ii++)
		{
			int i = solute_list[ii];
			float _R_tmp = pow( x[i-1][0] - coor_center[0] ,2) + \
						   pow( x[i-1][1] - coor_center[1] ,2) +\
						   pow(x[i-1][2] - coor_center[2], 2) ;
			if (_R_tmp > R_solute)
			{
				R_solute = _R_tmp;
			}
		}

		R_solute = sqrt(R_solute) + dmax;
		R_solute = pow(R_solute,2);

		// cout<<"Step 1 finished."<<endl;
		// #print "center %f\t%f\t%f" %(coor_x,coor_y,coor_z)
		// # print X_min,X_max
		// #step 2
		// #Get the range of the box.
		X_min = X_min - dmax;
		Y_min = Y_min - dmax;
		Z_min = Z_min - dmax;
		X_max = X_max + dmax;
		Y_max = Y_max + dmax;
		Z_max = Z_max + dmax;
		// # print X_min,X_max

		// # bin   = dmax *2.0 / nbins
		float bin = 0.05;
		int x_bins = floor((floor((X_max - X_min) /bin)/3)+1) *3 ; 
		int y_bins = floor((floor((Y_max - Y_min) /bin)/3)+1) *3 ;
		int z_bins = floor((floor((Z_max - Z_min) /bin)/3)+1) *3 ;

		// #step 4 
		// #assign each grid to solute atoms. 
		// #grid_in_solute contains that each grid blongs to which solute atom.


		map<int, pair<float,int> > grids_in_which_solute ; //=dict()
		map<int, vector<int> > solute_contain_grids ;//  =dict()

		for (int _x_temp=0; _x_temp < x_bins/3 ; _x_temp ++) // in range(x_bins/3)
		{
			// cout<< _x_temp<<endl;
			for (int _y_temp =0; _y_temp < y_bins/3; _y_temp ++)
			{
				// cout << _x_temp << _y_temp << endl;
				for (int _z_temp =0; _z_temp < z_bins/3; _z_temp ++)
				{

					vector<float> _grid_site;
					float grid_site[3];
					grid_site[0] =(X_min+(3*_x_temp+1.5)*bin);
					grid_site[1] =(Y_min+(3*_y_temp+1.5)*bin);
					grid_site[2] =(Z_min+(3*_z_temp+1.5)*bin);    


					float _test_1 = pow(coor_center[0]-grid_site[0],2)\
									+ pow(coor_center[1]-grid_site[1],2) \
									+ pow(coor_center[2]-grid_site[2],2);        
					if (_test_1 > R_solute)
					{
						continue;
					}    	
					// cout<< _grid_site
					vector<pair<float,int> > dist_index;
					dist_index.clear();
					// pair<float,int> dist_index
					dist_index = Sort_dist(x,solute_list,grid_site,coor_center);

					if (dist_index[0].first > dmax)
					{
						continue;
					}

					vector<int> _new_index;
					_new_index.clear();
					for (int i=0;i<15;i++)
					{
						_new_index.push_back(dist_index[i].second);
						// cout<<dist_index[i].first<<" ";
					}

					// cout<<endl;

					for (int _ii=0; _ii<3; _ii ++)
					{
						for (int _jj=0; _jj< 3; _jj ++)
						{
							for (int _kk=0; _kk <3; _kk ++)
							{
								// _grid_site.clear();
								grid_site[0]=(X_min+(3*_x_temp+_ii+0.5)*bin);
								grid_site[1]=(Y_min+(3*_y_temp+_jj+0.5)*bin);
								grid_site[2]=(Z_min+(3*_z_temp+_kk+0.5)*bin);

								float _test_1 = pow(coor_center[0]-grid_site[0],2)\
												+ pow(coor_center[1]-grid_site[1],2) \
												+ pow(coor_center[2]-grid_site[2],2);
								if (_test_1 > R_solute)
								{
									continue;
								}


								pair <float,int> min_dist_index; 
								min_dist_index= Min_dist(x,_new_index,grid_site,coor_center);
								// min_dist_index = dist_index[0];
								if (min_dist_index.second ==0)
								{
									continue;
								}

								int min_index = min_dist_index.second;

								int _serial = (3*_x_temp+_ii)*y_bins*z_bins + (3*_y_temp+_jj)*z_bins + (3*_z_temp+_kk);


								grids_in_which_solute[_serial]=min_dist_index;
								// cout << min_dist_index.first<< endl;

								if (solute_contain_grids.count(min_index) )
								{
									solute_contain_grids[min_index].push_back(_serial);
								}
								else
								{
									vector<int> new_vector;
									solute_contain_grids[min_index]=new_vector;
									solute_contain_grids[min_index].push_back(_serial);
								}
							}
						}
					}
				}
			}
		}


		// #step 5.
		// cout<<"step 5"<<endl;
		// #assign each solvent atom to grids.
		map<int, vector<int> > grid_in_solvent; //=[list() for i in range(x_bins * y_bins * z_bins)]
		for (int ii=0; ii< solvent_list.size(); ii++)
		{ 
			int i= solvent_list[ii]; //'' in solvent_list:
			float SV_x = x[i-1][0];
			float SV_y = x[i-1][1];
			float SV_z = x[i-1][2];
			int x_temp,  y_temp, z_temp;

			if (SV_x > X_min && SV_x < X_max)
			{
				x_temp = floor( (SV_x - X_min) / bin );
			}
			else
			{
				continue;
			}
			if (SV_y > Y_min && SV_y < Y_max)
			{    
				y_temp = floor( (SV_y - Y_min) / bin );
			}
			else
			{
				continue;
			}
			if (SV_z > Z_min && SV_z < Z_max)
			{
				z_temp = floor( (SV_z - Z_min) / bin );
			}
			else
			{
				continue;
			}
			int _serial = x_temp*y_bins*z_bins + y_temp*z_bins + z_temp; 
			// cout<< _serial<<endl;
			if (grid_in_solvent.count(_serial) )
			{
				grid_in_solvent[_serial].push_back(i);
			}
			else
			{
				vector<int> new_vector;
				grid_in_solvent[_serial]=new_vector;
				grid_in_solvent[_serial].push_back(i);

			}
		}
		// cout<< grid_in_solvent.size()<< endl;



		// # density   = MDAnalysis.core.units.convert(1.0, 'water', 'Angstrom^{-3}')*math.pow(10,3)
		// density   = MDAnalysis.core.units.convert(1.0, 'water', 'Angstrom^{-3}')
		// #        print density
		float density = 1.0 ; 
		float unit_conc =  (pow(bin,3) * density); //#unit solvent atom density.
		// cout<<unit_conc<< endl;

		vector<float> temp1  ; //  =list() //#A list used to contain grid_dist.
		vector<int>   temp2  ; //  =list() //#A list used to contain sol number for each grad.

		int TOTAL_ATOMS = 0;
		for (int ii =0;ii< solute_list.size();ii++)
		{
			int i = solute_list[ii];
			vector<int> temp_grids;
			if (solute_contain_grids.count( i) )
			{
				temp_grids=solute_contain_grids[i];
				// #print "solute %d, grids number %d" %(i,len(temp_grids))
			}
			else
			{
				continue;
			}
			// cout << temp_grids.size()<< endl;
			for (int ii=0;ii<temp_grids.size();ii++)// grid in temp_grids:
			{

				int grid = temp_grids[ii];
				// cout<< grid<< endl;
				int sol_number=grid_in_solvent[grid].size();
				pair<float,int> dist_index =grids_in_which_solute[grid];
				temp1.push_back(dist_index.first);
				// cout<< dist_index.first<< endl;
				temp2.push_back(sol_number);
			}
		}
		// cout << "hello world"<<endl;
		// cout << temp1.size()<< endl;

		vector<float> rdf_atom;
		rdf_atom=hist(temp1, temp2, 0.0, dmax, 100);
		// cout << "hello world"<< endl;
		for(int i=0;i<100;i++)
		{
			float add =0.0;
			add += rdf_atom[i];
			if (add > 0)
			{
				TOTAL_ATOMS += 1;
			}
		}
		// if sum(rdf_atom) > 0:
		//     TOTAL_ATOMS += 1
		for(int i=0;i<100;i++)
		{
			rdf_atom[i] = rdf_atom[i] / unit_conc;
		}
		// rdf_atom=numpy.array(rdf_atom) / unit_conc

		// # sys.exit()
		for(int i=0;i<100;i++)
		{
			GR[i] =GR[i] +rdf_atom[i];
		}
		// GR += rdf_atom;
		// #        print GR
		//         plt.clf()
		//         ax = plt.gca()
		//         ax.set_xlabel("Distance (nm)")
		//         ax.set_ylabel("pRDF(r)")
		//         x_label=[i*dmax/1000. for i in range(100)]
		// y_label=[GR[i]/ts.frame for i in range(100)]
		//         ax.plot(x_label,y_label,'b',)
		//         plt.draw()
		//         temp_filename="temp"+"%d"%(ts.frame)+".png"
		//         plt.savefig(temp_filename)
		for (int i=0; i<100;i++)
		{
			Data_file<<setw(8)<<setprecision(3)<<GR[i]/FRAMES;
		}
		Data_file<<"\n";
		Data_file.flush();
	}
	xdrfile_close(xtc);																																																	
}



// def Check_args():
// 	if len(sys.argv) != 5:
// 	print "Usage: pRDF_test.py coor_file traj_file index_file solute_index"
// 	else:
// 	coor_file    = sys.argv[1]
// 	traj_file    = sys.argv[2]
// 	index_file   = sys.argv[3]
// 	solute_index = int(sys.argv[4])
// 	pRDF(traj_file,coor_file,index_file,solute_index,10)

// 	if __name__ == '__main__':
// Check_args()


int main(int argc,char * argv[])
{
	// unit is nm
	char * coor_file;
	char * traj_file;
	char * index_file;
	char * data_file;


	switch(argc)
	{
		case 5:
			coor_file = argv[1];
			traj_file = argv[2];
			index_file = argv[3];
			data_file = argv[4];
			pRDF(traj_file,coor_file,index_file,data_file,1.0);
			break;

		case 2:
			Print_usage();
			exit(0);

		default:
			Print_usage();
			exit(0);
	}
}


void Print_usage()
{
	cout<<"\t\t"<<":-)  pRDF  (-:"<<endl;
	cout<<"------------------------------------------------"<<endl;
	cout<<"This is is a program designed for calculating proximal rdf of solvent."<<endl;
	cout<<endl;
	cout<<"------------------------------------------------"<<endl;
	cout << "Usage: pRDF coor_file traj_file index_file data_file"<< endl;
	cout<<"------------------------------------------------"<<endl;
	cout<<"Written by Zhu H. VERSION 1.0.0"<<endl;
	cout<<endl;
}

