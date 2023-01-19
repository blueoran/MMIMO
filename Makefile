OPT = 

.PHONY: debug all clean

all: mysim.exe 

mysim.exe: mysim.cpp decoder.cpp
	g++ $(OPT) -O3 -I/usr/include/eigen3 mysim.cpp -o mysim.exe

debug: mysim.cpp decoder.cpp
	g++ $(OPT) -O0 -g -I/usr/include/eigen3 mysim.cpp -o mysim-d.exe

clean:
	$(RM) *.exe