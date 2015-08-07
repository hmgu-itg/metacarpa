IDIR=/nfs/users/nfs_a/ag15/boost_1_55_0
CC=c++
LDIR=/nfs/users/nfs_a/ag15/boost_1_55_0/stage/lib
CFLAGS=-O3 -std=c++0x -static -I$(IDIR) -L$(LDIR)
LIBS=-lboost_program_options -lboost_serialization -lpthread

metacapa: metacapa3.cpp
	$(CC) $(CFLAGS) -o metacapa2 metacapa3.cpp $(LIBS)
