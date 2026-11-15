all: astral caster sister trial wtrial wcaster

dir:
	g++ -v 2>&1 | tail -n 1
	echo 'If installation failed, please ensure g++ version >= 13'
	mkdir -p bin

astral: dir
	g++ -std=c++20 -march=native -Ofast -D ASTRAL src/driver.cpp -o bin/astral

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

caster-dist: dir
	g++ -std=c++20 -march=native -Ofast src/caster-dist.cpp -o bin/caster-dist

clean: dir
	rm bin/*

doc: all
	mkdir -p doc
	bin/caster -H > doc/caster.md
	bin/astral -H > doc/astral.md
