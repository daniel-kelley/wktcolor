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

INC :=

WARN := -Wall
WARN += -Wextra
WARN += -Werror

CPPFLAGS := $(INC) -MP -MMD
CXXFLAGS := $(WARN) $(DEBUG)

SRC := wktcolor.cpp

OBJ := $(SRC:%.cpp=%.o)
DEP := $(SRC:%.cpp=%.d)

LDFLAGS :=

LDLIBS += -lboost_program_options
LDLIBS += -ligraph
LDLIBS += -lOpenMeshCore
LDLIBS += -lColPack
LDLIBS += -lgeos_c
LDLIBS += -lwkt

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

test: $(PROG) test_del.wkt
	./$(PROG) -gtest.off test_del.wkt > test.gml
	wktplot -TX -ctest.gml test_del.wkt

test_del.wkt: test_p.wkt
	wktdel $< > $@

test_p.wkt:
	wktrand  -q0.1 -r1.0 -x4 -y4 -n 10 > $@

clean:
	-rm -rf $(OBJ) $(DEP) $(PROG) \
		compile_commands.json $(CODECHECKER_OUT) \
		test*

-include $(DEP)
