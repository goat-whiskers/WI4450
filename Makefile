CXX=g++
CXX_FLAGS=-O2 -g -fopenmp -std=c++17

DEFS=-DNDEBUG

all: run_tests.x main_cg_poisson.x main_benchmarks.x

#default target (built when typing just "make")
#default: run_tests.x main_cg_poisson.x

# general rule to comple a C++ source file into an object file
%.o: %.cpp
	${CXX} -c ${CXX_FLAGS} ${DEFS} $<

#define some dependencies on headers
operations.o: operations.hpp timer.hpp
cg_solver.o: cg_solver.hpp operations.hpp timer.hpp
cg_poisson.o: cg_solver.hpp operations.hpp timer.hpp
gtest_mpi.o: gtest_mpi.hpp
main_benchmarks.o: operations.hpp timer.hpp 

TEST_SOURCES = test_operations.cpp test_cg_solver.cpp 
MAIN_OBJ = main_cg_poisson.o cg_solver.o operations.o timer.o
BENCH_OBJ = main_benchmarks.o operations.o timer.o

run_tests.x: run_tests.cpp ${TEST_SOURCES} gtest_mpi.o operations.o cg_solver.o
	${CXX} ${CXX_FLAGS} ${DEFS} -o run_tests.x $^

main_cg_poisson.x: main_cg_poisson.o ${MAIN_OBJ}
	${CXX} ${CXX_FLAGS} ${DEFS} -o main_cg_poisson.x $^

main_benchmarks.x: main_benchmarks.o ${BENCH_OBJ}
	${CXX} ${CXX_FLAGS} ${DEFS} -o main_benchmarks.x $^

test: run_tests.x
	./run_tests.x

bench: main_benchmarks.x
	./main_benchmarks.x

clean:
	-rm *.o *.x

# phony targets are run regardless of dependencies being up-to-date
PHONY: clean, test, bench
