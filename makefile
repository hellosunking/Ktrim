bin/ktrim: src/ktrim.cpp src/common.h src/util.h src/param_handler.h src/pe_handler.h src/se_handler.h
	cd src; g++ ktrim.cpp -fopenmp -O3 -o ../bin/ktrim; cd ..

install: bin/ktrim	# requires root
	cp bin/ktrim /usr/local/bin

clean:
	rm -f bin/ktrim

