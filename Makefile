CFLAGS = -lm

main: main.o
	gcc main.o -o main -lm

main.o: main.c
	gcc -c main.c -o main.o

run:
	./main

clean:
	rm -f main main.o