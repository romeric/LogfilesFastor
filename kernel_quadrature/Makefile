
# all: icc_ clang_ gcc_
all: clang_ gcc_


icc_:
	echo "Building and running FEM benchmarks with ICC"
	cd 2d_benchmark_mooney_rivlin;      ./run_benchmark.sh icc; cd ..
	cd 2d_benchmark_electro1;           ./run_benchmark.sh icc; cd ..
	cd 2d_benchmark_electro2;           ./run_benchmark.sh icc; cd ..

	echo "MooneyRivlin"
	cd 3d_benchmark_mooney_rivlin;      ./run_benchmark.sh icc; cd ..
	echo "Electromechanical Model 1"
	cd 3d_benchmark_electro1;           ./run_benchmark.sh icc; cd ..
	echo "Electromechanical Model 2"
	cd 3d_benchmark_electro2;           ./run_benchmark.sh icc; cd ..

clang_:
	echo "Building and running FEM benchmarks with CLANG"
	cd 2d_benchmark_mooney_rivlin;      ./run_benchmark.sh clang; cd ..
	cd 2d_benchmark_electro1;           ./run_benchmark.sh clang; cd ..
	cd 2d_benchmark_electro2;           ./run_benchmark.sh clang; cd ..

	echo "MooneyRivlin"
	cd 3d_benchmark_mooney_rivlin;      ./run_benchmark.sh clang; cd ..
	echo "Electromechanical Model 1"
	cd 3d_benchmark_electro1;           ./run_benchmark.sh clang; cd ..
	echo "Electromechanical Model 2"
	cd 3d_benchmark_electro2;           ./run_benchmark.sh clang; cd ..

gcc_:
	echo "Building and running FEM benchmarks with GCC"
	cd 2d_benchmark_mooney_rivlin;      ./run_benchmark.sh gcc; cd ..
	cd 2d_benchmark_electro1;           ./run_benchmark.sh gcc; cd ..
	cd 2d_benchmark_electro2;           ./run_benchmark.sh gcc; cd ..

	echo "MooneyRivlin"
	cd 3d_benchmark_mooney_rivlin;      ./run_benchmark.sh gcc; cd ..
	echo "Electromechanical Model 1"
	cd 3d_benchmark_electro1;           ./run_benchmark.sh gcc; cd ..
	echo "Electromechanical Model 2"
	cd 3d_benchmark_electro2;           ./run_benchmark.sh gcc; cd ..

clean:
	cd 2d_benchmark_mooney_rivlin;      make clean; cd ..
	cd 2d_benchmark_electro1;           make clean; cd ..
	cd 2d_benchmark_electro2;           make clean; cd ..

	cd 3d_benchmark_mooney_rivlin;      make clean; cd ..
	cd 3d_benchmark_electro1;           make clean; cd ..
	cd 3d_benchmark_electro2;           make clean; cd ..


.SILENT: