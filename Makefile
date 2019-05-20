CC=clang++
CFLAGS=-Wall -std=c++17 -stdlib=libc++

default:
	@mkdir -p bin
	make -C vendor/simdcomp
	$(CC) $(CFLAGS) -O3 -o bin/compress vendor/simdcomp/libsimdcomp.so.0.0.3 compress.cpp
	$(CC) $(CFLAGS) -O3 -o bin/decompress vendor/simdcomp/libsimdcomp.so.0.0.3 decompress.cpp

debug:
	@mkdir -p bin
	make -C vendor/simdcomp
	$(CC) $(CFLAGS) --debug -DDEBUG -o bin/compress vendor/simdcomp/libsimdcomp.so.0.0.3 compress.cpp
	$(CC) $(CFLAGS) --debug -DDEBUG -o bin/decompress vendor/simdcomp/libsimdcomp.so.0.0.3 decompress.cpp

.PHONY: clean
clean:
	@rm -r bin
