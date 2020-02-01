CC=clang++
# Can we take out zstd/lib/common or zstd/lib?
INCLUDES=-I. -Ivendor/zstd/lib -Ivendor/zstd/lib/common -Ivendor/zstd/lib/compress
CFLAGS=-Wall -std=c++17 -stdlib=libc++ $(INCLUDES)

ZSTD_FILES=$(wildcard vendor/zstd/lib/**/*.o)

default:
	@mkdir -p bin
	make -C vendor/simdcomp
	make lib-release -C vendor/zstd
	$(CC) $(CFLAGS) -O3 -o bin/compress vendor/simdcomp/libsimdcomp.so.0.0.3 $(ZSTD_FILES) compress.cpp
	$(CC) $(CFLAGS) -O3 -o bin/decompress vendor/simdcomp/libsimdcomp.so.0.0.3 $(ZSTD_FILES) decompress.cpp

debug:
	@mkdir -p bin
	make -C vendor/simdcomp
	make lib-release -C vendor/zstd
	$(CC) $(CFLAGS) -O0 --debug -DDEBUG -o bin/compress vendor/simdcomp/libsimdcomp.so.0.0.3 $(ZSTD_FILES) compress.cpp
	$(CC) $(CFLAGS) -O0 --debug -DDEBUG -o bin/decompress vendor/simdcomp/libsimdcomp.so.0.0.3 $(ZSTD_FILES) decompress.cpp

.PHONY: clean
clean:
	@rm -r bin
