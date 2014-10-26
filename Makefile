XDRINC = /home/zhuh/install/xdrfile/include/xdrfile/
XDRLIB = /home/zhuh/install/xdrfile/lib/
pRDF : main.o string_operate.o
	g++ -O2 -I${XDRINC} -o pRDF main.o string_operate.o -lxdrfile  -L${XDRLIB}
main.o : main.cpp pdb.h read_ndx.h string_operate.h 
	g++ -O2 -I${XDRINC} -c main.cpp
string_operate.o : string_operate.h string_operate.cpp
	g++ -O2 -c string_operate.cpp

clean :
	rm *.o
