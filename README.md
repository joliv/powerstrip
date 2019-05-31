# Powerstrip

```
make
bin/compress data.txt compressed.pstp
bin/decompress compressed.pstp data.txt
```

This works for files with one floating point value per line. Other modes are:

- `-b` binary: unsigned 16-bit binary integers
- `-i` text ints: one integer value per line

For example:

```
bin/compress -b data.bin compressed.pstp
bin/decompress -b compressed.pstp data.bin
```

Although the original and round-tripped formats do not need to match!
