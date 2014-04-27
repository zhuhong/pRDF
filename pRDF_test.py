#!/usr/bin/python
'''
TODO:
    write all the atoms of the solvent molecule to the index file.
'''
import math
import sys
import os
import numpy
import MDAnalysis
from itertools import izip
from MDPackage import Index
from MDPackage import Simple_atom
from MDPackage import usage
import MDPackage

import time as Time
import matplotlib.pyplot as plt


def Min_dist(solute_atoms,solu_l,solv_coor,coor_list,R_solute):
    '''
    atom_l: A full atom list.
    solu_l: A solute atom index list.
    solv_i: Index for one solvent molecule.
    '''

    _test_1 = (coor_list[0]-solv_coor[0])**2 \
    + (coor_list[1] - solv_coor[1])**2 \
    + (coor_list[2] - solv_coor[2])**2
    if _test_1 > R_solute:
        return 0,0

    dist_temp=[[abs(solute_atoms[i].atom_coor_x -solv_coor[0]),\
    abs(solute_atoms[i].atom_coor_y - solv_coor[1]),
    abs(solute_atoms[i].atom_coor_z - solv_coor[2])] for i in solu_l]

    # print dist_temp

    # min_dist  =min(dist_temp)
    # min_index =solu_l[dist_temp.index(min_dist)]
    min_dist = (coor_list[0]-solv_coor[0])**2 \
    + (coor_list[1] - solv_coor[1])**2 \
    + (coor_list[2] - solv_coor[2])**2
    min_dist = math.sqrt(min_dist)

    min_index =0
    # print min_dist
    # print solu_l
    for i in range(len(solu_l)):
        if dist_temp[i][0] < min_dist and dist_temp[i][1] < min_dist and dist_temp[i][2] < min_dist:
            tmp = math.sqrt((dist_temp[i][0])**2 + (dist_temp[i][1])**2 + (dist_temp[i][2])**2)
            if tmp < min_dist:
                min_dist = tmp
                min_index = i
                # print min_dist, min_index
            else:
                pass
        else:
            pass
    # if min_dist < 1:
        # print min_dist,min_index
    return min_dist,min_index

def Dist_Square(vect_1,vect_2):
    dist=(vect_2[0]-vect_1[0])**2 + (vect_2[1]-vect_1[1])**2 + (vect_2[2]-vect_1[2])**2
    return dist


def Get_solvent_list(atom_list):
    solvent_list = list()
    for atom in atom_list:
        if atom_list[atom].residue_name == "WAT" and atom_list[atom].atom_name == "O":
            solvent_list.append(atom)
        elif atom_list[atom].residue_name == "SOL" and atom_list[atom].atom_name == "OW":
            solvent_list.append(atom)
    return solvent_list

def hist(value_list,number_list,min_value, max_value, nbins=100):
    '''
    temp_list in the form [dist,number]
    '''
    # print value_list
    # print number_list
    # sys.exit()
    bin       = (max_value - min_value) / nbins
    
    grid_list =[0 for i in range(nbins)]
    numb_list =[0 for i in range(nbins)]

    for i in range(len(value_list)):
        temp = int( (value_list[i] - min_value) / bin )
        if value_list[i] < max_value:
            grid_list [temp] += 1
            numb_list[temp] += number_list[i]
        else:
       #     print "continued as vaule=%f" %value_list[i]
            continue

    #print grid_list
    #print numb_list
    for i in range(nbins):
        if grid_list[i] > 0:
            numb_list[i] /= float(grid_list[i])
        else:
            continue

    return numb_list


def pRDF(traj_file,coor_file,index_file,solute_index,dmax=20):
    '''
    A simple pRDF test here.
    '''

    START_TIME       =Time.time()

    HAS_TRAJFILE = False
    if os.path.isfile(traj_file):
        HAS_TRAJFILE = True
    
    atom_list    =Simple_atom.Get_atom_list(coor_file)
    index_list   =Index.Read_index_to_Inclass(index_file)
    solute_list  =index_list[solute_index].group_list

    solvent_list =Get_solvent_list(atom_list)

    if HAS_TRAJFILE:
        U=MDAnalysis.Universe(coor_file,traj_file)
    else:
        U=MDAnalysis.Universe(coor_file)

    GR=numpy.zeros((100),dtype=numpy.float64)
    EXIST_NUM = 0
    if os.path.isfile("datafile.xvg"):
        Data_file = open("datafile.xvg",'r+')
        lines = Data_file.readlines()
        EXIST_NUM = len(lines)
        try:
            last_line = lines[-1].split()
            for x in range(100):
                GR[x] = float(last_line[x])*EXIST_NUM
        except:
            pass
    else:
        Data_file=open("datafile.xvg",'w')
    #step 1 
    # Get the center of the solute. 
    for ts in U.trajectory:

        print "Checking frame number: %d" %(ts.frame)
        if ts.frame < EXIST_NUM +1:
            continue
        solute_atoms = dict()
        
        coor_x=0.0
        coor_y=0.0
        coor_z=0.0
        X_min = ts._x[0]
        X_max = ts._x[0]
        Y_min = ts._y[0]
        Y_max = ts._y[0]
        Z_min = ts._z[0]
        Z_max = ts._z[0]        
        R_solute =0.0

        for i in solute_list:
            coor_x += ts._x[i-1]
            coor_y += ts._y[i-1]
            coor_z += ts._z[i-1]
            u_atom = MDPackage.Coor.unit_atom.unit_atom(atom_coor_x=ts._x[i-1],atom_coor_y=ts._y[i-1],atom_coor_z=ts._z[i-1])
            solute_atoms[i]=u_atom

            if ts._x[i-1] < X_min:
                X_min = ts._x[i-1]
            elif ts._x[i-1] > X_max:
                X_max = ts._x[i-1]

            if ts._y[i-1] < Y_min:
                Y_min = ts._y[i-1]
            elif ts._y[i-1] > Y_max:
                Y_max = ts._y[i-1]

            if ts._z[i-1] < Z_min:
                Z_min = ts._z[i-1]
            elif ts._z[i-1] > Z_max:
                Z_max = ts._z[i-1]

        coor_x /= len(solute_list)
        coor_y /= len(solute_list)
        coor_z /= len(solute_list)

        for i in solute_list:
            _R_tmp = ( solute_atoms[i].atom_coor_x - coor_x ) **2 + \
            ( solute_atoms[i].atom_coor_y - coor_y ) **2 +\
            ( solute_atoms[i].atom_coor_z - coor_z ) **2 
            if _R_tmp > R_solute:
                R_solute = _R_tmp

        R_solute = math.sqrt(R_solute) + dmax
        R_solute = R_solute **2
        # print R_solute
        # sys.exit()


        # print "Step 1 finished."
        #print "center %f\t%f\t%f" %(coor_x,coor_y,coor_z)
        # print X_min,X_max
        #step 2
        #Get the range of the box.
        X_min = X_min - dmax
        Y_min = Y_min - dmax
        Z_min = Z_min - dmax
        X_max = X_max + dmax
        Y_max = Y_max + dmax
        Z_max = Z_max + dmax
        # print X_min,X_max

        # bin   = dmax *2.0 / nbins
        bin = 1.0
        x_bins = int((X_max - X_min) /bin) +1
        y_bins = int((Y_max - Y_min) /bin) +1
        z_bins = int((Z_max - Z_min) /bin) +1 

        #print "bin:",bin

        #step 4 
        #assign each grid to solute atoms. 
        #grid_in_solute contains that each grid blongs to which solute atom.


        grids_in_which_solute =dict()
        solute_contain_grids  =dict()

        # print x_bins,y_bins,z_bins
        _gauss_value = -1
        for i in range(x_bins * y_bins * z_bins ):
            z_temp = i / ( x_bins * y_bins) 
            y_temp = (i % (x_bins * y_bins)) /  x_bins 
            x_temp = i % x_bins

            grid_site=[X_min+(x_temp+0.5)*bin, Y_min+(y_temp+0.5)*bin, Z_min+(z_temp+0.5)*bin]
            #print grid_site
            min_dist, min_index= Min_dist(solute_atoms,solute_list,grid_site,[coor_x,coor_y,coor_z],R_solute)
            if min_index == 0:
                continue
            # _gauss_value = min_dist
            if i % 300 ==0:
                NOW_TIME=Time.time()
                BIN_TIME=NOW_TIME-START_TIME
                sys.stderr.write("grid ID %10d, time used: %6.2f s\r" %(i,BIN_TIME))
                sys.stderr.flush()
            try:
                grids_in_which_solute[i]=[min_index,min_dist]
            except:
                print "hello to see you"
            try:
                solute_contain_grids[min_index].append(i)
            except:
                solute_contain_grids[min_index]=list()
                solute_contain_grids[min_index].append(i)
                #print solute_contain_grids

        # for i in solute_contain_grids:
        #     print i,len(solute_contain_grids[i])
        # sys.exit()

        # print "\nStep 4 finished."

        #step 5.
        #assign each solvent atom to grids.
        grid_in_solvent=[list() for i in range(x_bins * y_bins * z_bins)]
        for i in solvent_list:
            SV_x = ts._x[i-1]
            SV_y = ts._y[i-1]
            SV_z = ts._z[i-1]

            if SV_x > X_min and SV_x < X_max:
                x_temp = int( (SV_x - X_min) / bin )
            else:
                continue
            if SV_y > Y_min and SV_y < Y_max:    
                y_temp = int( (SV_y - Y_min) / bin )
            else:
                continue
            if SV_z > Z_min and SV_z < Z_max:
                z_temp = int( (SV_z - Z_min) / bin )
            else:
                continue
            grid_in_solvent[z_temp*x_bins*y_bins + y_temp*x_bins + x_temp].append(i)
            # print "append solvent %d" %i

        # print "Step 5 finished."
         # step 6.
        #calculating the g(r) for each solute atom.

       # density   = MDAnalysis.core.units.convert(1.0, 'water', 'Angstrom^{-3}')*math.pow(10,3)
        density   = MDAnalysis.core.units.convert(1.0, 'water', 'Angstrom^{-3}')
#        print density
        unit_conc =  ((bin)**3 * density) #unit solvent atom density.
#        print unit_conc

        temp1     =list() #A list used to contain grid_dist.
        temp2     =list() #A list used to contain sol number for each grad.
 
        # TOTAL_ATOMS = 0
        for i in solute_list:
            try:
                temp_grids=solute_contain_grids[i]
                #print "solute %d, grids number %d" %(i,len(temp_grids))
            except:
                continue
        #rdf_atom=[0 for i in range(50)]
            # bin       =float(dmax)/nbins


        
            for grid in temp_grids:

                sol_number=len(grid_in_solvent[grid])
                # if sol_number ==0:
                #     continue
                #print "   %d" %sol_number,
                # try:
                blabla,dist=grids_in_which_solute[grid]
                # except:
                    # continue
            
                temp1.append(dist)
                temp2.append(sol_number)

            
            # if len(temp1) == 0:
            #     continue
                # print temp1
                # print temp2
        #if i == 10:
        #    sys.exit()
        rdf_atom=hist(temp1, temp2, 0.0, dmax, 100)
        # print rdf_atom
        # print unit_conc
        #print rdf_atom
        # if sum(rdf_atom) > 0:
        #     TOTAL_ATOMS += 1
        rdf_atom=numpy.array(rdf_atom) / unit_conc

        # sys.exit()
        GR += rdf_atom
#        print GR
        plt.clf()
        ax = plt.gca()
        ax.set_xlabel("Distance (nm)")
        ax.set_ylabel("pRDF(r)")
        x_label=[i*dmax/1000. for i in range(100)]
        y_label=[GR[i]/ts.frame for i in range(100)]
        ax.plot(x_label,y_label,'b',)
        plt.draw()
        temp_filename="temp"+"%d"%(ts.frame)+".png"
        plt.savefig(temp_filename)
        for i in range(100):
            Data_file.write("%12.4f" %y_label[i])
        Data_file.write("\n")
        Data_file.flush()
        

    GR = GR / U.trajectory.numframes
    Data_file.close()
#   print TOTAL_ATOMS
#        print len(solute_index)
    for i in range(32):
        print "%f\t%f" %(2.0/32*(i+0.5),GR[i])

  #  print GR


def Check_args():
    if len(sys.argv) != 5:
        print "Usage: pRDF_test.py coor_file traj_file index_file solute_index"
    else:
        coor_file    = sys.argv[1]
        traj_file    = sys.argv[2]
        index_file   = sys.argv[3]
        solute_index = int(sys.argv[4])
        pRDF(traj_file,coor_file,index_file,solute_index,10)

if __name__ == '__main__':
    plt.ion()
    Check_args()

