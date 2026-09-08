# Physics regression tests.  They compile against the OpenACC headers, which are
# plain host C++, so no GPU or CUDA toolkit is needed:
#
#     make -C tests run

CXX ?= g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wno-unknown-pragmas -I../OpenACC_cli

TESTS := geodesic_invariants metric_codegen_equiv

all: $(TESTS)

%: %.cpp ../OpenACC_cli/cuda_ray.h ../OpenACC_cli/black_hole.h
	$(CXX) $(CXXFLAGS) -o $@ $< -lm

run: $(TESTS)
	@for t in $(TESTS); do echo "== $$t"; ./$$t || exit 1; done

clean:
	$(RM) $(TESTS)

.PHONY: all run clean
