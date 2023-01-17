all: mydecoder.so

mydecoder.so: Decoder.cpp ../decoder.cpp
	g++ --shared $< -O3 -fPIC -I/usr/include/eigen3 -o $@

clean:
	$(RM) *.so

.PHONY: all clean