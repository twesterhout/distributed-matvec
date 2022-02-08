.RULES:
.POSIX:

OPTIMIZATION ?= --debug
CFLAGS = -Ithird_party/include $(OPTIMIZATION)
HDF5_LIBS ?= `pkg-config --libs hdf5` -lhdf5_hl
LDFLAGS += -Lthird_party/lib -llattice_symmetries_haskell -llattice_symmetries $(HDF5_LIBS) -lutil -lgomp -lpthread

all: basis

# plugin.o: plugin.c plugin.h
# libplugin.a: plugin.o
# 	$(AR) rcs $@ $^

test: test.chpl
	chpl $(OPTIMIZATION) -o $@ $^

# basis: basis.chpl states.chpl
# 	chpl $(CFLAGS) --main-module basis -o $@ $^ $(LDFLAGS) 

# basis: lib/libbasis.a
# lib/libbasis.a: basis.chpl states.chpl
# 	chpl $(CFLAGS) --library --static --library-makefile -o basis $^ $(LDFLAGS) 
DATA_FILES = data/heisenberg_chain_10.h5 \
	     data/heisenberg_kagome_12.h5 \
	     data/heisenberg_kagome_16.h5 \
	     data/heisenberg_square_4x4.h5

.PHONY: check_io
check_io: test_basis_construction test_vector_loading $(DATA_FILES)
	for numLocales in 4 3 2 1; do \
	  for yamlPath in \
	    "heisenberg_chain_10.yaml" \
	    "heisenberg_kagome_12.yaml" \
	    "heisenberg_kagome_16.yaml" \
	    "heisenberg_square_4x4.yaml"; do \
	    ./test_vector_loading -nl $${numLocales} \
	      --kInputDataPath "data/$${yamlPath%.yaml}.h5" \
	      --kOutputDataPath "output.h5"; \
	    ./test_basis_construction -nl $${numLocales} \
	      --kInputBasisPath "data/$${yamlPath}" \
	      --kInputDataPath "data/$${yamlPath%.yaml}.h5" \
	      --kOutputDataPath "output.h5"; \
	  done; \
	done

test_matvec: test_matvec.chpl matvec.chpl Distribute.chpl basis.chpl states.chpl io.chpl merge.chpl wrapper.chpl
	CHPL_TARGET_CPU=native chpl \
		$(CFLAGS) \
		-o $@ $^ \
		--main-module $@ \
		$(LDFLAGS)

test_basis_construction: test_basis_construction.chpl Distribute.chpl basis.chpl states.chpl io.chpl merge.chpl wrapper.chpl
	CHPL_TARGET_CPU=native chpl \
		$(CFLAGS) \
		-o $@ $^ \
		--main-module $@ \
		$(LDFLAGS)

test_vector_loading: test_vector_loading.chpl basis.chpl states.chpl io.chpl merge.chpl wrapper.chpl
	CHPL_TARGET_CPU=native chpl \
		$(CFLAGS) \
		-o $@ $^ \
		--main-module $@ \
		$(LDFLAGS)

data/%.h5: data/%.yaml data/SpinED
	cd data && \
	OMP_NUM_THREADS=`nproc` ./SpinED $(<F)

data/SpinED:
	mkdir -p data && cd data && \
	curl -LJO https://github.com/twesterhout/spin-ed/releases/download/manual/SpinED-4c3305a && \
	mv SpinED-4c3305a SpinED && \
	chmod +x SpinED

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
