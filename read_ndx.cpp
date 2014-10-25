#include "read_ndx.h"
#include "string_operate.h"



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
        printf("Group %4d (%15s) has %6d elements\n", i , index_list[i].group_name , index_list[i].group_list.size())
    }
}
// vector<string>  read_2_group(char * index_file)
// {
// 	vector<string> vs;
// 	string s;
// 	ifstream in(index_file);
// 	while(!in.fail())
// 	{
// 		getline(in,s);
// 		if(s.find("[")<10)
// 		{
// 			vs.push_back(Split(Split(Split(s,'[',0),']',0),' ',0));
// 		}
// 	}
// 	in.close();
// 	return vs;
// }

// vector<int>  read_2_group_n(char * index_file)
// {
// 	vector<int> vs;
// 	string s;
// 	int num=0;
// 	ifstream in(index_file);
// 	getline(in,s);
// 	while(!in.fail())
// 	{	
// 		while(s.find("[")<10)
// 		{
// 			getline(in,s);
// 			while(s.find("[")>100&&s.size()>0)
// 			{
// 				num+=Split_v(s,' ').size();
// 				getline(in,s);
// 			}
// 			vs.push_back(num);
// 			num=0;
// 		}
// 		getline(in,s);
// 	}
// 	in.close();
// 	return vs;
// }


// vector<int>  read_group_2_serial(char * index_file, string group_name)
// {
// 	vector<int> vi;
// 	string s;
// 	ifstream in(index_file);
// 	while(!in.fail())
// 	{
// 		getline(in,s);
// 		if((s.find("[")<10)&(Split(Split(Split(s,'[',0),']',0),' ',0)==group_name))
// 		{
// 			getline(in,s);
// 			vector<string> temp;
// 			while(s.find("[")>10&&!in.fail())
// 			{
// 				temp.clear();
// 				temp=Split_v(s,' ');
// 				for(int i=0;i<temp.size();i++)
// 				{
// 					vi.push_back(atoi(temp.at(i).c_str()));
// 				}
// 				getline(in,s);
// 			}
			
// 		}
// 	}
// 	in.close();
// 	return vi;
// }