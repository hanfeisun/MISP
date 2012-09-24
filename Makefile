SOURCES = motif.c mis.c
CC = gcc
LIB = -lz -lm


all: kseq.h mis.c motif.c
	$(CC) -O3 -g $(SOURCES) -o mis $(LIB)

dbg: kseq.h mis.c motif.c
	$(CC) -g $(SOURCES) -o mis $(LIB)

clean:
	rm -f *.o mis
	tm test_all

check-syntax:
	$(CC) -Wall -Wextra -pedantic -fsyntax-only $(SOURCES) $(LIB)


test:
	make all
	./mis test.seq database/cistrome.db 0.001 all test
