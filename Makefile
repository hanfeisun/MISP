SOURCES = motif.c mis.c
CC = gcc

all: kseq.h mis.c motif.c
		$(CC) -g -O2 $(SOURCES) -o mis -lz

dbg: kseq.h mis.c motif.c
		$(CC) -g $(SOURCES) -o mis -lz 
clean:
		rm -f *.o mis

check-syntax:
		$(CC) -Wall -Wextra -pedantic -fsyntax-only $(SOURCES)
