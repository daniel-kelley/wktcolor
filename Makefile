#
#  Makefile
#
#  Copyright (c) 2021 by Daniel Kelley
#

DEBUG ?= -g

PREFIX ?= /usr/local

# address thread undefined etc.
ifneq ($(SANITIZE),)
DEBUG += -fsanitize=$(SANITIZE)
endif

INC ?=
CXXFLAGS := $(INC) -MP -MMD

WARN := -Wall
WARN += -Wextra
WARN += -Werror
CPPFLAGS := $(WARN) $(DEBUG)

SRC := wktcolor.cpp

OBJ := $(SRC:%.cpp=%.o)
DEP := $(SRC:%.cpp=%.d)

LDFLAGS := $(DEBUG)
LDFLAGS += -lboost_program_options
# for CGAL?
#LDFLAGS += -lgmp
#LDFLAGS += -lmpfr

LDFLAGS += -lOpenMeshCore
LDFLAGS += -lgeos_c
LDFLAGS += -lwkt

#
# CodeChecker support
#
BEAR ?= bear
CODECHECKER ?= CodeChecker
CODECHECKER_OUT = ./reports
CODECHECKER_ARGS := --ctu
CODECHECKER_ARGS += --skip .codechecker.skip
CODECHECKER_ARGS += --output $(CODECHECKER_OUT)
#CODECHECKER_ARGS += --verbose debug

include .codechecker.mk

PROG := wktcolor

.PHONY: all install uninstall check clean

all: $(PROG)

$(PROG): $(SRC)

install: $(PROG)
	install -p -m 755 $(PROG) $(PREFIX)/bin

uninstall:
	-rm -f $(PREFIX)/bin/$(PROG)

	-rm -f $(PROG) $(OBJ) $(DEP)

check: compile_commands.json
	-$(CODECHECKER) analyze $< $(CODECHECKER_ARGS)
	$(CODECHECKER) parse $< $(CODECHECKER_OUT)

compile_commands.json:
	$(BEAR) $(MAKE)

clean:
	-rm -rf $(OBJ) $(DEP) $(PROG) \
		compile_commands.json $(CODECHECKER_OUT)


-include $(DEP)
