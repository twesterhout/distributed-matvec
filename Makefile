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

test_basis_construction: test_basis_construction.chpl basis.chpl states.chpl io.chpl merge.chpl wrapper.chpl
	# CHPL_TARGET_CPU=native chpl $(CFLAGS) --fast --vectorize -o $@ $^ $(LDFLAGS) 
	CHPL_TARGET_CPU=native chpl \
		-Ithird_party/include \
		--debug -o $@ $^ \
		--main-module $@ \
		-Lthird_party/lib \
		-llattice_symmetries_haskell \
		-llattice_symmetries \
		`pkg-config --libs hdf5` -lhdf5_hl \
		-lutil \
		-lgomp \
		-lpthread

basis: basis.chpl states.chpl merge.chpl
	# CHPL_TARGET_CPU=native chpl $(CFLAGS) --fast --vectorize -o $@ $^ $(LDFLAGS) 
	CHPL_TARGET_CPU=native chpl \
		-Ithird_party/include \
		--fast -o $@ $^ \
		-Lthird_party/lib \
		-llattice_symmetries_haskell \
		-llattice_symmetries \
		`pkg-config --libs hdf5` -lhdf5_hl \
		-lutil \
		-lgomp \
		-lpthread

copying: copying.chpl
	CHPL_TARGET_CPU=native chpl --fast -s algorithm=1 -o $@ $^

merge: merge.chpl
	chpl $(CFLAGS) --debug -o merge $^ $(LDFLAGS) 

junk: junk.chpl
	chpl --debug -o $@ $^

main: main.c
	gcc -o $@ \
		-Ithird_party/include \
		$< \
		-pthread \
		-fopenmp \
		-Lthird_party/lib \
		-llattice_symmetries_haskell \
		-llattice_symmetries \
		-lgmp \
		-lnuma \
		-lm \
		-lz \
		-ldl \
		-lutil

# -include lib/Makefile.basis
# main: main.c basis
# 	$(CHPL_COMPILER) $(CHPL_CFLAGS) -o $@ $< $(CHPL_LDFLAGS)

# Dummy file we use to reproduce internal compiler errors in Chapel for
# submitting issues.
.PHONY: error
error: error.chpl
	chpl -o $@ $^

.PHONY: dependencies
dependencies: assets/third_party.tar.bz2
	tar -xf $<

ghc-8.10.7.sif: ghc-8.10.7.def
	singularity build --fakeroot --force --no-cleanup $@ $<

01_ghc.sif: 01_ghc.def
	singularity build --fakeroot --force $@ $<

02_hdf5.sif: 02_hdf5.def
	singularity build --fakeroot --force $@ $<

assets/hdf5.tar.bz2: 02_hdf5.sif
	singularity run --bind assets:/prefix $<

03_lattice-symmetries.sif: 03_lattice-symmetries.def
	singularity build --fakeroot --force $@ $<

assets/lattice-symmetries.tar.bz2: 03_lattice-symmetries.sif
	singularity run --bind assets:/prefix $<

04_haskell.sif: 04_haskell.def assets/hdf5.tar.bz2 assets/lattice-symmetries.tar.bz2
	singularity build --fakeroot --force $@ $<

05_dependencies.sif: 05_dependencies.def 04_haskell.sif
	singularity build --fakeroot --force --no-cleanup $@ $<

assets/third_party.tar.bz2: 05_dependencies.sif
	singularity run --bind assets:/prefix $<

.PHONY: clean
clean:
	rm -rf lib/ basis_server_real main
