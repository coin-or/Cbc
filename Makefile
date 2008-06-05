COININSTDIR = $(HOME)/COIN/Osi/0.98/build
CXX = g++
CXXFLAGS = -g -Wall -I$(COININSTDIR)/include/coin

NAME = dynamicbranching
SRC	= dynamicbranching.cpp
OBJS 	= $(SRC:.cpp=.o)

LDFLAGS  = -Wl,-rpath,-L$(COININSTDIR)/lib
LDFLAGS += -L$(COININSTDIR)/lib
LDFLAGS += -lOsiClp 
LDFLAGS += -lOsi 
LDFLAGS += -lClp 
LDFLAGS += -lCoinUtils 
LDFLAGS += `cat $(COININSTDIR)/share/doc/coin/Osi/osi_addlibs.txt`
LDFLAGS += `cat $(COININSTDIR)/share/doc/coin/Clp/clp_addlibs.txt`
LDFLAGS += `cat $(COININSTDIR)/share/doc/coin/CoinUtils/coinutils_addlibs.txt`


$(NAME): $(OBJS) $(SRC) $(INCL)  
	$(CXX) $(OBJS) $(LDFLAGS) -o $(NAME) 

clean:
	rm -f $(NAME) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< 
