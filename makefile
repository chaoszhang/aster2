all: caster

dir:
	g++ -v 2>&1 | tail -n 1
	echo 'If installation failed, please ensure g++ version >= 13'
	mkdir -p bin
	
caster: dir
	g++ -std=c++20 -march=native -Ofast -D CASTER src/driver.cpp -o bin/caster

doc: all
	mkdir -p doc
	bin/caster -H > doc/caster.md