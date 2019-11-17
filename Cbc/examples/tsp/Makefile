#
# Makefile for TSP Example in COIN-OR CBC
# haroldo (unxusr), 2019-11-11 16:46

CXXFLAGS=-Og -g -c `pkg-config --cflags cbc`
LDFLAGS=-Og -g `pkg-config --libs cbc`

all:tsp-subtour
	

tsp-subtour:tsp-subtour.o
	g++ $(LDFLAGS) tsp-subtour.o -o tsp-subtour

tsp-subtour.o:tsp-subtour.cpp
	g++ $(CXXFLAGS) tsp-subtour.cpp -o tsp-subtour.o

clean:
	rm -f tsp-subtour *.o


# vim:ft=make
#
