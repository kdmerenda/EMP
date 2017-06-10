CC=g++ -g -fopenmp -O2
OBJS=emp2d.o
SRC=emp2d-fork.cpp
all: emp2d
emp2d: $(OBJS)
	$(CC) -o emp2d $(OBJS)
emp2d.o: $(SRC)
	$(CC) -o emp2d.o -c $(SRC)
clean:
	rm *.o emp2d
