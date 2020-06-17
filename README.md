<p align="center"><img src="powerstrip.png" align="center" width="300px"></p>

Powerstrip is an integer compression algorithm designed for energy meter data, particularly for device-level or circuit-level data. It provides high compression ratios, rapid compression and decompression speeds, and nearly zero lossiness.

See the [Powerstrip paper](https://www.bowdoin.edu/~sbarker/research/pdf/eenergy20-powerstrip.pdf) for details, which appeared at [ACM e-Energy 2020](https://energy.acm.org/conferences/eenergy/2020/).

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
