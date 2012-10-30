
SOURCES = motif.c mis.c
CC = cc
LIB = -lz -lm


all: kseq.h mis.c motif.c 
	$(CC) -Wall -g $(SOURCES) -o mis $(LIB)

dbg: kseq.h mis.c motif.c
	$(CC) -g $(SOURCES) -o mis $(LIB)

clean:
	rm -f *.o mis
	rm test_*
	rm test2_*

check-syntax:
	$(CC) -Wall -Wextra -pedantic -fsyntax-only $(SOURCES) $(LIB)


test:
	make all
	time ./mis test.seq database/cistrome.db 0.001 all test
	time ./mis test.seq database/cistrome.db 0.001 EN0055 test
	time ./mis test2.seq database/cistrome.db 0.001 all test2
