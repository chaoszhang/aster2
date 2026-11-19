all: astral astral-iad caster sister trial wtrial wcaster qdist4p

dir:
	g++ -v 2>&1 | tail -n 1
	echo 'If installation failed, please ensure g++ version >= 13'
	mkdir -p bin

astral: dir
	g++ -std=c++20 -march=native -Ofast -D ASTRAL src/driver.cpp -o bin/astral

astral-iad: dir
	g++ -std=c++20 -march=native -Ofast -D ASTRAL_IAD src/driver.cpp -o bin/astral-iad

caster: dir
	g++ -std=c++20 -march=native -Ofast -D CASTER src/driver.cpp -o bin/caster

sister: dir
	g++ -std=c++20 -march=native -Ofast -D SISTER src/driver.cpp -o bin/sister

trial: dir
	g++ -std=c++20 -march=native -Ofast -D TRIAL src/driver.cpp -o bin/trial

wtrial: dir
	g++ -std=c++20 -march=native -Ofast -D WTRIAL src/driver.cpp -o bin/wtrial

wcaster: dir
	g++ -std=c++20 -march=native -Ofast -D WCASTER src/driver.cpp -o bin/wcaster

qdist4p: dir
	g++ -std=c++20 -march=native -Ofast -D QDIST4P src/driver.cpp -o bin/qdist4p

caster-dist: dir
	g++ -std=c++20 -march=native -Ofast src/caster-dist.cpp -o bin/caster-dist

clean: dir
	rm bin/*

doc: all
	mkdir -p doc
	bin/caster -H > doc/caster.md
	bin/astral -H > doc/astral.md
