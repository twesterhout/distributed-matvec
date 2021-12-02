.POSIX:

CFLAGS = -I/home/tom/src/lattice-symmetries-haskell/cbits `pkg-config --cflags lattice_symmetries`
LDFLAGS = -L/home/tom/src/lattice-symmetries-haskell/build -llattice_symmetries_haskell `pkg-config --libs lattice_symmetries` 

all: basis

# plugin.o: plugin.c plugin.h
# libplugin.a: plugin.o
# 	$(AR) rcs $@ $^

# test: test.chpl states.chpl libplugin.a
# 	chpl -o $@ -L. -lplugin $(LDFLAGS) plugin.h $<

# basis: basis.chpl states.chpl
# 	chpl $(CFLAGS) --main-module basis -o $@ $^ $(LDFLAGS) 

# basis: lib/libbasis.a
# lib/libbasis.a: basis.chpl states.chpl
# 	chpl $(CFLAGS) --library --static --library-makefile -o basis $^ $(LDFLAGS) 

basis: basis.chpl states.chpl merge.chpl
	chpl $(CFLAGS) --debug -o basis $^ $(LDFLAGS) 

merge: merge.chpl
	chpl $(CFLAGS) --debug -o merge $^ $(LDFLAGS) 

junk: junk.chpl
	chpl --debug -o $@ $^

# -include lib/Makefile.basis
# main: main.c basis
# 	$(CHPL_COMPILER) $(CHPL_CFLAGS) -o $@ $< $(CHPL_LDFLAGS)

# Dummy file we use to reproduce internal compiler errors in Chapel for
# submitting issues.
.PHONY: error
error: error.chpl
	chpl -o $@ $^

# lib/libstates.a: states.chpl libplugin.a
# 	chpl --static --library --library-makefile -L. -lplugin $(LDFLAGS) plugin.h $<

# foo: foo.chpl
# 	chpl --static --library --library-makefile $<

# include lib/Makefile.foo

# bar: bar.c lib/libfoo.a
# 	$(CHPL_COMPILER) $(CHPL_CFLAGS) -o $@ $< $(CHPL_LDFLAGS)

.PHONY: clean
clean:
	rm -rf lib/ basis_server_real main
