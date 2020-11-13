objects = xdrfile.o xdrfile_xtc.o main.o FCorrelatorSAt.o
a.out : $(objects)
	g++ -o a.out $(objects)
xdrfile.o : xdrfile.h
	gcc -c xdrfile.c
xdrfile_xtc.o : xdrfile_xtc.h xdrfile.h
	gcc -c xdrfile_xtc.c
FCorrelatorSAt.o : FCorrelatorSAt.h
	g++ -c FCorrelatorSAt.cc
main.o : xdrfile.h xdrfile_xtc.h FCorrelatorSAt.h
	g++ -c main.cc
.PHONY : clean
clean :
	rm a.out $(objects)