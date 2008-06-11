COININSTDIR = $(HOME)/COIN/Osi/0.98/build
CXX = g++
CXXFLAGS = -g -Wall -I$(COININSTDIR)/include/coin

LDFLAGS  = -Wl,-rpath,$(COININSTDIR)/lib
LDFLAGS += -L$(COININSTDIR)/lib
LDFLAGS += -lOsiClp 
LDFLAGS += -lOsi 
LDFLAGS += -lClp 
LDFLAGS += -lCoinUtils 
LDFLAGS += `cat $(COININSTDIR)/share/doc/coin/Osi/osi_addlibs.txt`
LDFLAGS += `cat $(COININSTDIR)/share/doc/coin/Clp/clp_addlibs.txt`
LDFLAGS += `cat $(COININSTDIR)/share/doc/coin/CoinUtils/coinutils_addlibs.txt`


dynamicbranching: dynamicbranching.o
	$(CXX) dynamicbranching.o $(LDFLAGS) -o $@

dynamicLL: dynamicLL.o
	$(CXX) dynamicbranching.o $(LDFLAGS) -o $@

clean:
	rm -f *.o dynamicbranching dynamicLL

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< 
