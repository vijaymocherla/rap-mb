cd rap_mb/
icpx -O3 -fPIC -Wall -shared -std=c++17 py_rap_mb.cc rap_mb.cc -I/usr/local/include -L/usr/local/lib \
       	-lsundials_arkode -lsundials_nvecserial -lsundials_core \
	-qopenmp -lstdc++ -qmkl \
	-undefined dynamic_lookup $(python3 -m pybind11 --includes) \
	-o _rap_mb$(python3-config --extension-suffix)

icpx -O3 -std=c++17 test.cc rap_mb.cc -I/usr/local/include -L/usr/local/lib \
       	-lsundials_arkode -lsundials_nvecserial -lsundials_core \
		-lm -qopenmp -lstdc++ -qmkl \
		-o test.x
