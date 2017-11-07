CC := g++
CCFLAGS := -g -std=c++11

all : euler velvert carparr

euler : code-Euler.cpp
	$(CC) $(CCFLAGS) $< -o $@

velvert : code-VelVert.cpp
	$(CC) $(CCFLAGS) $< -o $@

carparr : code-CarParr.cpp
	$(CC) $(CCFLAGS) $< -o $@

.phony : clean

clean:
	rm *.o euler velvert


