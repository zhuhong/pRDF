#ifndef READ_NDX_H
#define READ_NDX_H

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>



#include "string_operate.h"

using namespace std;

class Index_class
{
    public:
     	string group_name;
     	vector<int> group_list;
};


vector<Index_class> Read_index_to_Inclass(char * filename)
{
        ifstream fp(filename);
        vector<Index_class> vs;
        Index_class temp_Index_class;
        string s;

    while(!fp.fail())
    {
    	getline(fp,s);
        if (s.find("[")<10)
        {
        	if(temp_Index_class.group_name.length()>0)
        	{
        		Index_class temp = temp_Index_class;
        		vs.push_back(temp);
        		Index_class temp_Index_class;
        		
        	}
        	temp_Index_class.group_name= Split(Split(Split(s,'[',0),']',0),' ',0);
        	temp_Index_class.group_list.clear();
        }
        else if (s.length()>1)
        {
        	vector<string> temp;
			temp.clear();
			temp=Split_v(s,' ');
			for(int i=0;i<temp.size();i++)
			{
				temp_Index_class.group_list.push_back(atoi(temp.at(i).c_str()));

			}
			// cout<< temp_Index_class.group_list.size()<<endl;
       	}
    }
        vs.push_back(temp_Index_class);
        fp.close();
        return vs;
}


void Print_Index(vector<Index_class> index_list)
{
    for(int i=0;i< index_list.size();i++)
    {
        printf("Group %4d (%15s) has %6d elements\n", i , index_list[i].group_name.c_str() , index_list[i].group_list.size());
    }
}


#endif