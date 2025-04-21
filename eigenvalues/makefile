OPTS = -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
a.out: main.o functions.o solve.o
	g++ $(OPTS) main.o functions.o solve.o -o a.out
main.o: main.cpp functions.h
	g++ -c $(OPTS) main.cpp -lm
functions.o: functions.cpp functions.h
	g++ -c $(OPTS) functions.cpp
solve.o: solve.cpp functions.h
	g++ -c $(OPTS) solve.cpp
clean:
	rm -f main.o functions.o solve.o a.out
