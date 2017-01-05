

INC = -I../../../../ -I../
OPT = -std=c++11 -O3 -mavx

.SILENT:


all: $(bench)


fastor:
	$(CXX) benchmark_quadrature.cpp $(INC) $(OPT) $(PDEG) $(CHECKS)
	./a.out $(FOLDER)

scalar_variants:
	$(CXX) benchmark_quadrature_scalar.cpp $(INC) $(OPT) $(PDEG) $(TCROSS) $(VECT) $(CHECKS)
	./a.out $(FOLDER)

clean:
	rm -rf a.out out_* results_*