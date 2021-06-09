# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

GRCHOMBO_SOURCE = $(shell pwd)/../../Source

ebase := InterpolatorTest

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/GRChomboCore \
			$(GRCHOMBO_SOURCE)/utils \
			$(GRCHOMBO_SOURCE)/simd \
			$(GRCHOMBO_SOURCE)/AMRInterpolator

include $(CHOMBO_HOME)/mk/Make.test
