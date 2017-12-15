exefiles=hot-icc-O1 hot-icc-O2 hot-icc-O2-avx hot-icc-g-O2 hot-icc-g
all:hot-icc-O1 hot-icc-O2 hot-icc-O2-avx hot-icc-g-O2 hot-icc-g
hot-icc-O1:hot-intel.cpp
	icc -O1 -I./3rd-lib-base/gsl2.4/include -I./Eigen3 -L./3rd-lib-base/gsl2.4/lib hot-intel.cpp -o hot-icc-O1 -lgsl -lgslcblas
hot-icc-O2:hot-intel.cpp
	icc -O2 -I./3rd-lib-base/gsl2.4/include -I./Eigen3 -L./3rd-lib-base/gsl2.4/lib hot-intel.cpp -o hot-icc-O2 -lgsl -lgslcblas
hot-icc-O2-avx:hot-intel.cpp
	icc -O2 -xavx -qopt-report=5 -I./3rd-lib-base/gsl2.4/include -I./Eigen3 -L./3rd-lib-base/gsl2.4/lib hot-intel.cpp -o hot-icc-O2-avx -lgsl -lgslcblas
hot-icc-g-O2:hot-intel.cpp
	icc -g -O2 -I./3rd-lib-base/gsl2.4/include -I./Eigen3 -L./3rd-lib-base/gsl2.4/lib hot-intel.cpp -o hot-icc-g-O2 -lgsl -lgslcblas
hot-icc-g:hot-intel.cpp
	icc -g -fp-model precise -I./3rd-lib-base/gsl2.4/include -I./Eigen3 -L./3rd-lib-base/gsl2.4/lib hot-intel.cpp -o hot-icc-g -lgsl -lgslcblas
.phony:clean
clean:
	rm -rf ${exefiles}
