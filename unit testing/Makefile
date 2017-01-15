CC = gcc
OBJECTS = main.o rmsd.o drmsd.o clrec.o nnlsh.o
CFLAGS = -lm -lblas -llapacke

recommendation: $(OBJECTS)
	$(CC) -o recommendation $(OBJECTS) $(CFLAGS)

main.o: main.c
	$(CC) -c main.c 

rmsd.o: rmsd.c
	$(CC) -c rmsd.c 

drmsd.o: drmsd.c
	$(CC) -c drmsd.c
	
clrec.o: clrec.c
	$(CC) -c clrec.c

nnlsh.o: nnlsh.c
	$(CC) -c nnlsh.c

.PHONY: clean

clean:
	rm -f recommendation $(OBJECTS)
