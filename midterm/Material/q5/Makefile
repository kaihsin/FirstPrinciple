CC := g++
CCFLAGS := -g -std=c++11

all : euler velvert predcorr

euler : code-Euler.cpp
	$(CC) $(CCFLAGS) $< -o $@

velvert : code-VelVert.cpp
	$(CC) $(CCFLAGS) $< -o $@

predcorr : code-PredCorr.cpp
	$(CC) $(CCFLAGS) $< -o $@

.phony : clean

clean:
	rm *.o euler velvert predcorr


