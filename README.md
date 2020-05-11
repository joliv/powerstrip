# Powerstrip ⚡

Quick-start for an out-of-source build:

```bash
git clone --recursive https://github.com/joliv/powerstrip
mkdir build # out-of-source build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

This results in the two binaries: `build/pstrip` and `build/unpstrip`, as well as the library `build/libpowerstrip.*`.

Currently, we only compress binary files of unsigned 16-bit integers:

```bash
# Compress input.bin to output.pstrip
pstrip input.bin output.pstrip

# Decompress output.pstrip back to roundtripped.bin
unpstrip output.pstrip roundtripped.bin
```
