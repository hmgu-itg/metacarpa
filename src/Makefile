#IDIR=/nfs/users/nfs_a/ag15/boost_1_55_0
#IDIR=/nfs/team144/software/boost_1.6/boost_1_60_0
IDIR=/usr/include
CC=c++
#LDIR=/nfs/users/nfs_a/ag15/boost_1_55_0/stage/lib
#LDIR=/nfs/team144/software/boost_1.6/boost_1_60_0/stage/lib
LDIR=lib
CFLAGS=-O3 -std=c++11 -I$(IDIR) -L$(LDIR)
LIBS=-lboost_program_options -lboost_serialization -lpthread

metacapa: metacarpa.cpp
	$(CC) $(CFLAGS) -o metacarpa metacarpa.cpp $(LIBS)
