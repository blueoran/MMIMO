.PHONY: all clean

all: mysim.exe 

mysim.exe: mysim.cpp decoder.cpp
	g++ -O3 -I/usr/include/eigen3 mysim.cpp -o mysim.exe

clean:
	$(RM) mysim.exe