# Makefile for bench/fj_bench
#
# Builds against the installed MIPster library via pkg-config.
# Run from this directory:
#   make
#   make PREFIX=~/prog/cbc   (or set PREFIX in environment)

PREFIX ?= $(HOME)/prog/cbc
PKG    := $(PREFIX)/lib/pkgconfig/mipster.pc

CXX      ?= g++
CXXFLAGS := -O2 -std=c++17 $(shell PKG_CONFIG_PATH=$(PREFIX)/lib/pkgconfig pkg-config --cflags mipster)
# feasibilityjump.hh is in the installed include dir
INCSRC   := -I$(PREFIX)/include/mipster
LDFLAGS  := $(shell PKG_CONFIG_PATH=$(PREFIX)/lib/pkgconfig pkg-config --libs mipster)

# Add the monorepo src dir so we can include feasibilityjump.hh directly
# (it's not installed as a public header since it's an internal third-party file).
# Fall back gracefully if building outside the source tree.
SRCDIR   := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))../src
ifneq ($(wildcard $(SRCDIR)/feasibilityjump.hh),)
  INCSRC += -I$(SRCDIR)
endif

CXXFLAGS += $(INCSRC)

TARGET := fj_bench

.PHONY: all clean

all: $(TARGET)

$(TARGET): fj_bench.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Built $@"

clean:
	rm -f $(TARGET)
