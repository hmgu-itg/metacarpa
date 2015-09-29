# Makefile assumes BOOST is installed as a system library in /usr/lib 
# and /usr/include. 
# If you have BOOST installed elsewhere, set BOOST on the make command line
# to the prefix in which BOOST is installed 
# (i.e. `make BOOST=/home/alice/boost_1_55_0`)
BOOST=/usr

# If you have not installed boost but have a built source tree, you will need 
# to instead set BOOST_INCLUDE and BOOST_LIB
# (i.e. `make BOOST_INCLUDE=/home/alice/boost_1_55_0_src BOOST_LIB=/home/alice/boost_1_55_0_src/stage/lib`)
BOOST_INCLUDE=$(BOOST)/include
BOOST_LIB=$(BOOST)/lib

# For static binary, set SHARED=0 (i.e. `make SHARED=0`)
SHARED=1

CC=c++
RM=rm

BOOST_LIBS=-L$(BOOST_LIB) -lboost_program_options -lboost_serialization
BOOST_CFLAGS=-I$(BOOST_INCLUDE)
CFLAGS=-O3 -std=c++11 
LIBS=-lpthread

ifeq ($(SHARED),0)
# For static binary, pass -static
	CFLAGS+=-static
else
# For dynamic binary, add rpath to BOOST_LIB directory
	BOOST_LIBS+=-Wl,-rpath=$(BOOST_LIB)
endif

metacarpa: metacarpa.cpp
	$(CC) $(CFLAGS) $(BOOST_CFLAGS) -o metacarpa metacarpa.cpp $(BOOST_LIBS) $(LIBS)

clean:
	$(RM) -f metacarpa

