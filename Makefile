.RULES:
.POSIX:

UNAME = $(shell uname)
OPTIMIZATION ?= --debug
COMPILER ?= llvm
CFLAGS = -Ithird_party/include $(OPTIMIZATION) --target-compiler=$(COMPILER)
LDFLAGS += -Lthird_party/lib -llattice_symmetries_haskell -llattice_symmetries_core
# ifeq ($(UNAME), Linux)
#   LDFLAGS += -lnuma
# endif

PRIMME_CFLAGS = -I/home/tom/src/primme/include
PRIMME_LDFLAGS = -L/home/tom/src/primme/lib -lprimme -llapacke -lopenblas -lm -lgomp -lpthread
HDF5_CFLAGS = $(shell pkg-config --cflags hdf5)
HDF5_LIBS = -lhdf5_hl $(shell pkg-config --libs hdf5)

ifeq ($(UNAME), Linux)
  CONDA_CC ?= $(shell conda run -n ci_devel bash -c "which \$${CC}")
  CONDA_PREFIX ?= $(shell conda run -n ci_devel bash -c "echo \$${CONDA_PREFIX}")
  SHARED_EXT = so
  SHARED_FLAG = -shared -rdynamic
else
  CONDA_CC = $(CC)
  CONDA_PREFIX = 
  SHARED_EXT = dylib
  SHARED_FLAG = -dynamiclib
endif

PREFIX = $(PWD)
PACKAGE = lattice-symmetries-chapel
GIT_COMMIT = $(shell git rev-parse --short HEAD)
DIST = $(PACKAGE)-$(GIT_COMMIT)

# MODULES = src/ApplyOperator.chpl src/StatesEnumeration.chpl src/helper.c
LIB_MODULES = src/LatticeSymmetries.chpl \
              src/FFI.chpl \
              src/ForeignTypes.chpl \
              src/Vector.chpl \
              src/ConcurrentAccessor.chpl \
              src/StatesEnumeration.chpl \
              src/BatchedOperator.chpl \
              src/DistributedMatrixVector.chpl \
              src/CommonParameters.chpl

APP_MODULES = $(LIB_MODULES) \
	      src/HashedToBlock.chpl \
	      src/BlockToHashed.chpl \
	      src/MyHDF5.chpl

.PHONY: all
all: examples

ifeq ($(UNAME), Linux)
  CHPL_LIBS = LD_LIBRARY_PATH=$(PWD)/third_party/lib:$$LD_LIBRARY_PATH
endif
ifeq ($(UNAME), Darwin)
  CHPL_LIBS = DYLD_LIBRARY_PATH=$(PWD)/third_party/lib:$$DYLD_LIBRARY_PATH
endif
CHPL_ARGS =

.PHONY: examples
examples: bin/Example02 bin/Example05

.PHONY: test
test: bin/TestStatesEnumeration bin/TestMatrixVectorProduct

.PHONY: check
check: check-states-enumeration check-matrix-vector-product

.PHONY: benchmark-states-enumeration
benchmark-states-enumeration: bin/TestStatesEnumeration
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_pyrochlore_2x2x2.yaml --kRepresentatives data/large-scale/construction/heisenberg_pyrochlore_2x2x2.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_36.yaml --kRepresentatives data/large-scale/construction/heisenberg_kagome_36.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_square_6x6.yaml --kRepresentatives data/large-scale/construction/heisenberg_square_6x6.h5

.PHONY: check-states-enumeration
check-states-enumeration: bin/TestStatesEnumeration
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_4.yaml --kRepresentatives data/matvec/heisenberg_chain_4.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_6.yaml --kRepresentatives data/matvec/heisenberg_chain_6.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_8.yaml --kRepresentatives data/matvec/heisenberg_chain_8.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_10.yaml --kRepresentatives data/matvec/heisenberg_chain_10.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_12.yaml --kRepresentatives data/matvec/heisenberg_chain_12.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_16.yaml --kRepresentatives data/matvec/heisenberg_chain_16.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_20.yaml --kRepresentatives data/matvec/heisenberg_chain_20.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_24.yaml --kRepresentatives data/matvec/heisenberg_chain_24.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_24_symm.yaml --kRepresentatives data/matvec/heisenberg_chain_24_symm.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_12.yaml --kRepresentatives data/matvec/heisenberg_kagome_12.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_12_symm.yaml --kRepresentatives data/matvec/heisenberg_kagome_12_symm.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_16.yaml --kRepresentatives data/matvec/heisenberg_kagome_16.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_4x4.yaml --kRepresentatives data/matvec/heisenberg_square_4x4.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_5x5.yaml --kRepresentatives data/matvec/heisenberg_square_5x5.h5

.PHONY: benchmark-matrix-vector-product
benchmark-matrix-vector-product: bin/TestMatrixVectorProduct
	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_pyrochlore_2x2x2.yaml --kRepresentatives data/heisenberg_pyrochlore_2x2x2.h5
	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_36.yaml --kRepresentatives data/heisenberg_kagome_36.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_6x6.yaml --kVectors data/large-scale/matvec/heisenberg_square_6x6.h5

.PHONY: check-matrix-vector-product
check-matrix-vector-product: bin/TestMatrixVectorProduct
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_4.yaml --kVectors data/matvec/heisenberg_chain_4.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_6.yaml --kVectors data/matvec/heisenberg_chain_6.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_8.yaml --kVectors data/matvec/heisenberg_chain_8.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_10.yaml --kVectors data/matvec/heisenberg_chain_10.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_12.yaml --kVectors data/matvec/heisenberg_chain_12.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_16.yaml --kVectors data/matvec/heisenberg_chain_16.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_20.yaml --kVectors data/matvec/heisenberg_chain_20.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_24.yaml --kVectors data/matvec/heisenberg_chain_24.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_24_symm.yaml --kVectors data/matvec/heisenberg_chain_24_symm.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_12.yaml --kVectors data/matvec/heisenberg_kagome_12.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_12_symm.yaml --kVectors data/matvec/heisenberg_kagome_12_symm.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_16.yaml --kVectors data/matvec/heisenberg_kagome_16.h5


TEST_DATA_URL = https://surfdrive.surf.nl/files/index.php/s/OK5527Awfgl1hT2/download

.PHONY: data/construction
data/construction:
	mkdir -p data && cd data && \
	wget -q -O tmp.zip $(TEST_DATA_URL)?path=%2Fdata%2Fconstruction && \
	unzip tmp.zip && rm tmp.zip

.PHONY: data/matvec
data/matvec:
	mkdir -p data && cd data && \
	wget -q -O tmp.zip $(TEST_DATA_URL)?path=%2Fdata%2Fmatvec && \
	unzip tmp.zip && rm tmp.zip

.PHONY: data/large-scale
data/large-scale:
	mkdir -p data && cd data && \
	wget -q -O tmp.zip $(TEST_DATA_URL)?path=%2Fdata%2Flarge-scale && \
	unzip tmp.zip && rm tmp.zip

lib: lib/liblattice_symmetries_chapel.so

lib/liblattice_symmetries_chapel.a: $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) --library --static --library-makefile -o lattice_symmetries_chapel $^ $(LDFLAGS)

lib/liblattice_symmetries_chapel.so: lib/liblattice_symmetries_chapel.a
	$(CONDA_CC) $(SHARED_FLAG) -o lib/liblattice_symmetries_chapel.$(SHARED_EXT) src/library.c $^ `$$CHPL_HOME/util/config/compileline --libraries` $(LDFLAGS)
# ifeq ($(UNAME), Darwin)
# 	install_name_tool -id lib/liblattice_symmetries_chapel.$(SHARED_EXT) lib/liblattice_symmetries_core.$(SHARED_EXT)
# endif

.PHONY: release
release: lib
	mkdir -p $(DIST)/include
	mkdir -p $(DIST)/lib
	install -m644 -C third_party/include/*.h $(DIST)/include/
	install -m644 -C third_party/lib/*.$(SHARED_EXT) $(DIST)/lib/
	install -m644 -C lib/liblattice_symmetries_chapel.$(SHARED_EXT) $(DIST)/lib/
ifeq ($(UNAME), Linux)
	find $(DIST)/lib/ -name "*.$(SHARED_EXT)" -exec patchelf --set-rpath '$$ORIGIN' {} \;
endif
	tar -cf $(DIST).tar $(DIST)
	rm -f $(DIST).tar.bz2
	bzip2 $(DIST).tar
ifneq ($(realpath $(PREFIX)), $(PWD))
	install -m644 -C $(DIST).tar.bz2 $(PREFIX)
endif
	rm -r $(DIST)

bin/TestStatesEnumeration: test/TestStatesEnumeration.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LIBS) $(LDFLAGS)

bin/TestMatrixVectorProduct: test/TestMatrixVectorProduct.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LIBS) $(LDFLAGS)

bin/Example01: example/Example01.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LIBS) $(LDFLAGS)

bin/Example02: example/Example02.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LIBS) $(LDFLAGS)

bin/Example03: example/Example03.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LIBS) $(LDFLAGS)

bin/Example04: example/Example04.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LIBS) $(LDFLAGS)

bin/Example05: example/Example05.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LIBS) $(LDFLAGS) 

bin/dummy: src/dummy.chpl
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ $<

bin/primme: src/PRIMME.chpl
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(PRIMME_CFLAGS) -o $@ $< $(PRIMME_LDFLAGS)

bin/Diagonalize: src/Diagonalize.chpl src/PRIMME.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(PRIMME_CFLAGS) -o $@ --main-module $(@F) $^ $(LDFLAGS) $(PRIMME_LDFLAGS)

# Dummy file we use to reproduce internal compiler errors in Chapel for
# submitting issues.
.PHONY: error
error: src/error.chpl
	chpl -o $@ $^

.PHONY: clean
clean:
	rm -rf bin/
