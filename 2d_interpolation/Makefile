CFLAGS = -O3 -pthread -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

a.out:	main.o functions.o solution.o residual.o
		g++ $^ $(CFLAGS) -o $@

main.o: main.cpp header.h
		g++ -c $(CFLAGS) main.cpp

residual.o: residual.cpp header.h
		g++ -c $(CFLAGS) residual.cpp

solution.o: solution.cpp header.h
		g++ -c $(CFLAGS) solution.cpp

functions.o: functions.cpp header.h
		g++ -c $(CFLAGS) functions.cpp

clean:
		rm -f *.out *.o