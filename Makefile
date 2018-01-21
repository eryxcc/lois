CFLAGS?=-std=c++0x -O3
CC=g++

all: out/tutorial.txt out/autotest.txt out/mintest.txt 

bin/tutorial: include/lois.h tests/tutorial.cpp bin/liblois.a
	$(CC) $(CFLAGS) tests/tutorial.cpp bin/liblois.a -o bin/tutorial

bin/liblois.a: obj/element.o obj/env.o obj/formula.o obj/relations.o obj/solvers.o obj/utils.o obj/extra.o obj/pseudo.o
	ar rvs bin/liblois.a obj/env.o obj/formula.o obj/relations.o obj/solvers.o obj/utils.o obj/element.o obj/extra.o obj/pseudo.o

obj/env.o: include/lois.h include/loisinternal.h src/env.cpp
	$(CC) $(CFLAGS) src/env.cpp -c -o obj/env.o

obj/element.o: include/lois.h include/loisinternal.h src/element.cpp
	$(CC) $(CFLAGS) src/element.cpp -c -o obj/element.o

obj/extra.o: include/lois.h include/loisextra.h src/extra.cpp
	$(CC) $(CFLAGS) src/extra.cpp -c -o obj/extra.o

obj/utils.o: include/lois.h include/loisinternal.h src/utils.cpp
	$(CC) $(CFLAGS) src/utils.cpp -c -o obj/utils.o

obj/formula.o: include/lois.h include/loisinternal.h src/formula.cpp
	$(CC) $(CFLAGS) src/formula.cpp -c -o obj/formula.o

obj/relations.o: include/lois.h include/loisinternal.h src/relations.cpp
	$(CC) $(CFLAGS) src/relations.cpp -c -o obj/relations.o

obj/solvers.o: include/lois.h include/loisinternal.h src/solvers.cpp
	$(CC) $(CFLAGS) src/solvers.cpp -c -o obj/solvers.o

obj/pseudo.o: include/lois.h include/loisextra.h src/pseudo.cpp
	$(CC) $(CFLAGS) src/pseudo.cpp -c -o obj/pseudo.o

bin/autotest: tests/autotest.cpp bin/liblois.a
	$(CC) $(CFLAGS) tests/autotest.cpp bin/liblois.a -o bin/autotest

out/autotest.txt: bin/autotest
	bin/autotest > out/autotest.txt

view-autotest: bin/autotest
	bin/autotest

out/soltest.txt: bin/soltest
	bin/soltest > out/soltest.txt

bin/mintest: tests/mintest.cpp bin/liblois.a
	$(CC) $(CFLAGS) tests/mintest.cpp -o bin/mintest bin/liblois.a
	
view-mintest: bin/mintest
	bin/mintest

bin/orbtest: tests/orbtest.cpp bin/liblois.a
	$(CC) $(CFLAGS) tests/orbtest.cpp -o bin/orbtest bin/liblois.a


bin/pushdown: tests/pushdown.cpp bin/liblois.a
	$(CC) $(CFLAGS) tests/pushdown.cpp -o bin/pushdown bin/liblois.a

bin/soltest: tests/soltest.cpp bin/liblois.a
	$(CC) $(CFLAGS) tests/soltest.cpp -o bin/soltest bin/liblois.a

bin/stackdemo: bin/liblois.a tests/stackdemo.cpp
	$(CC) $(CFLAGS) tests/stackdemo.cpp -o bin/stackdemo bin/liblois.a

out/mintest.txt: bin/mintest
	bin/mintest > out/mintest.txt

testout: bin/autotest
	bin/autotest

out/tutorial.txt: bin/tutorial
	bin/tutorial > out/tutorial.txt

bin/xtest: include/lois.h tests/xtest.cpp bin/liblois.a
	$(CC) $(CFLAGS) tests/xtest.cpp bin/liblois.a -o bin/xtest

view-xtest: bin/xtest
	bin/xtest

view-tutorial: bin/tutorial
	bin/tutorial

bin/learning: tests/learning.cpp bin/liblois.a
	$(CC) $(CFLAGS) tests/learning.cpp bin/liblois.a -o bin/learning

bin/csp: tests/csp.cpp bin/liblois.a
	$(CC) $(CFLAGS) tests/csp.cpp bin/liblois.a -o bin/csp

out/learning.txt: bin/learning
	bin/learning > out/learning.txt

bin/parsetutorial: manual/parsetutorial.cpp
	g++ -std=c++0x -O3 manual/parsetutorial.cpp -o bin/parsetutorial

manual/snippets/ok: bin/parsetutorial out/tutorial.txt
	bin/parsetutorial out/tutorial.txt tests/tutorial.cpp

manual/manual.pdf: manual/manual.tex manual/tech-intro.tex manual/syntax.tex manual/safety.tex manual/relations.tex manual/piecewise.tex manual/overview.tex manual/contents.tex manual/lois.bib manual/snippets/ok manual/tutorial.tex
	cd manual; pdflatex manual; bibtex manual; pdflatex manual; pdflatex manual


export: lois.tgz

lois.tgz: Makefile out/tutorial.txt out/autotest.txt tests/soltest.cpp 
	rm -rf lois/*
	mkdir -p lois lois/include lois/tests lois/src lois/obj lois/bin lois/out
	cp include/*.h lois/include/
	cp src/*.cpp lois/src/
	cp tests/tutorial.cpp tests/autotest.cpp tests/mintest.cpp tests/stackdemo.cpp tests/soltest.cpp lois/tests/
	cp Makefile lois/Makefile
	cp LICENSE lois/LICENSE
	tar zcf lois.tgz lois/ --owner=1000 --group=1000 --numeric-owner

lois-results.tgz: lois.tgz
	cd lois; make; cd ..
	tar zcf lois-results.tgz lois/ --owner=1000 --group=1000 --numeric-owner
