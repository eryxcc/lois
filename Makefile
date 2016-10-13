all: out/tutorial.txt out/autotest.txt out/mintest.txt

bin/tutorial: include/lois.h tests/tutorial.cpp bin/liblois.a
	g++ -std=c++0x -O3 tests/tutorial.cpp bin/liblois.a -o bin/tutorial

bin/liblois.a: obj/env.o obj/formula.o obj/relations.o obj/solvers.o obj/utils.o obj/element.o obj/extra.o obj/pseudo.o
	ar rvs bin/liblois.a obj/env.o obj/formula.o obj/relations.o obj/solvers.o obj/utils.o obj/element.o obj/extra.o obj/pseudo.o

obj/env.o: include/lois.h include/loisinternal.h src/env.cpp
	g++ -std=c++0x -O3 src/env.cpp -c -o obj/env.o

obj/element.o: include/lois.h include/loisinternal.h src/element.cpp
	g++ -std=c++0x -O3 src/element.cpp -c -o obj/element.o

obj/extra.o: include/lois.h include/loisextra.h src/extra.cpp
	g++ -std=c++0x -O3 src/extra.cpp -c -o obj/extra.o

obj/utils.o: include/lois.h include/loisinternal.h src/utils.cpp
	g++ -std=c++0x -O3 src/utils.cpp -c -o obj/utils.o

obj/formula.o: include/lois.h include/loisinternal.h src/formula.cpp
	g++ -std=c++0x -O3 src/formula.cpp -c -o obj/formula.o

obj/relations.o: include/lois.h include/loisinternal.h src/relations.cpp
	g++ -std=c++0x -O3 src/relations.cpp -c -o obj/relations.o

obj/solvers.o: include/lois.h include/loisinternal.h src/solvers.cpp
	g++ -std=c++0x -O3 src/solvers.cpp -c -o obj/solvers.o

obj/pseudo.o: include/lois.h include/loisextra.h src/pseudo.cpp
	g++ -std=c++0x -O3 src/pseudo.cpp -c -o obj/pseudo.o

bin/autotest: tests/autotest.cpp bin/liblois.a
	g++ -std=c++0x -O3 tests/autotest.cpp bin/liblois.a -o bin/autotest

out/autotest.txt: bin/autotest
	bin/autotest > out/autotest.txt

out/soltest.txt: bin/soltest
	bin/soltest > out/soltest.txt

bin/mintest: tests/mintest.cpp bin/liblois.a
	g++ -std=c++0x -O3 tests/mintest.cpp -o bin/mintest bin/liblois.a
	
bin/orbtest: tests/orbtest.cpp bin/liblois.a
	g++ -std=c++0x -O3 tests/orbtest.cpp -o bin/orbtest bin/liblois.a


bin/pushdown: tests/pushdown.cpp bin/liblois.a
	g++ -std=c++0x -O3 tests/pushdown.cpp -o bin/pushdown bin/liblois.a


bin/soltest: tests/soltest.cpp bin/liblois.a
	g++ -std=c++0x -O3 tests/soltest.cpp -o bin/soltest bin/liblois.a

bin/stackdemo: bin/liblois.a tests/stackdemo.cpp
	g++ -std=c++0x -O3 tests/stackdemo.cpp -o bin/stackdemo bin/liblois.a

out/mintest.txt: bin/mintest
	bin/mintest > out/mintest.txt

testout: bin/autotest
	bin/autotest

out/tutorial.txt: bin/tutorial
	bin/tutorial > out/tutorial.txt

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
