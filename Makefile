all: axle-demo


%.o: %.c
	gcc -c $^ -O3


axle-demo: axle.o doubles.o
	gcc -o $@ $^ -lm


clean:
	rm -f *.o axle-demo *~
