WORKDIR = 'pwd'
CC = g++
CFLAGS = -c
INC = -Iinc

gdpm : main.o tree.o io.o 
	$(CC) -o gdpm main.o tree.o io.o 

main.o : src/main.cpp
	$(CC) -c src/main.cpp $(INC)
tree.o : src/tree.cpp 
	$(CC) -c src/tree.cpp $(INC)
io.o : src/io.cpp 
	$(CC) -c src/io.cpp $(INC)

clean :
	rm gdpm \
	main.o tree.o io.o
