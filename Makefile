.RULES:
.POSIX:

UNAME = $(shell uname)
OPTIMIZATION ?= --debug
COMPILER ?= llvm
CFLAGS = -Ithird_party/include $(OPTIMIZATION) --target-compiler=$(COMPILER)
LDFLAGS += -Lthird_party/lib -llattice_symmetries_haskell -llattice_symmetries_core
ifeq ($(UNAME), Linux)
  LDFLAGS += -lnuma
endif

PRIMME_CFLAGS = -I/home/tom/src/primme/include
PRIMME_LDFLAGS = -L/home/tom/src/primme/lib -lprimme -lopenblas -lm -lgomp -lpthread

MODULES = src/ApplyOperator.chpl src/StatesEnumeration.chpl src/helper.c


bin/Example01: example/Example01.chpl $(MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ --main-module $(@F) $^ $(LDFLAGS)


# Dummy file we use to reproduce internal compiler errors in Chapel for
# submitting issues.
.PHONY: error
error: src/error.chpl
	chpl -o $@ $^


.PHONY: clean
clean:
	rm -rf bin/
