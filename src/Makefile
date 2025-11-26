# Detect OS and set correct shared library extension
UNAME := $(shell uname -s)

ifeq ($(UNAME), Darwin)  # macOS
  EXT = so
else ifeq ($(UNAME), Linux)
  EXT = so
else  # Assume Windows (MSYS or Cygwin)
  EXT = dll
endif

# Define all target shared object files
TMB_TARGETS = Gaussian_ind.$(EXT) Gaussian_ind_fixed.$(EXT) Gaussian_dep.$(EXT) Gaussian_dep_fixed.$(EXT) Poisson_ind.$(EXT) Poisson_ind_fixed.$(EXT)

# Default target: compile all shared objects
all: $(TMB_TARGETS)

# General rule to compile any .cpp file into a shared object
%.$(EXT): %.cpp
	$(R_HOME)/bin$(R_ARCH_BIN)/Rscript --vanilla -e "TMB::compile('$<', '', safebounds = FALSE,safeunload = FALSE)"

# Clean up generated shared object and intermediate files
clean:
	rm -f $(TMB_TARGETS) *.o
