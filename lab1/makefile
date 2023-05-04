SOURCE = g++ -Wall -lm -O2 -std=c++1z 

%.o: %.cpp
	$(SOURCE) -c $< -o $@

lab1.1: 1.1.o ./lib/matrix.o ./lib/functions.o
	$(SOURCE) 1.1.o ./lib/matrix.o ./lib/functions.o -o lab1.1

lab1.2: 1.2.o ./lib/matrix.o ./lib/functions.o
	$(SOURCE) 1.2.o ./lib/matrix.o ./lib/functions.o -o lab1.2

clean:
	rm *.o ./lib/*.o lab1.1 lab1.2
