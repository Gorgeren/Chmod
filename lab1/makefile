SOURCE = g++ -Wall -lm -O2 -std=c++1z -Wextra -g

%.o: %.cpp
	$(SOURCE) -c $< -o $@

lab1.1: 1.1.o ./lib/matrix.o ./lib/functions.o
	$(SOURCE) 1.1.o ./lib/matrix.o ./lib/functions.o -o lab1.1

lab1.2: 1.2.o ./lib/matrix.o ./lib/functions.o
	$(SOURCE) 1.2.o ./lib/matrix.o ./lib/functions.o -o lab1.2

lab1.3: 1.3.o ./lib/matrix.o ./lib/functions.o
	$(SOURCE) 1.3.o ./lib/matrix.o ./lib/functions.o -o lab1.3

lab1.4: 1.4.o ./lib/matrix.o ./lib/functions.o
	$(SOURCE) 1.4.o ./lib/matrix.o ./lib/functions.o -o lab1.4
clean:
	rm *.o ./lib/*.o lab1.1 lab1.2 lab1.3 lab1.4
