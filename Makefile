all: main.o matrix.o
	g++ -o m *.o
main.o:
	g++ -c main.cpp
matrix.o:
	g++ -c matrix.cpp
clean:
	rm *.o m