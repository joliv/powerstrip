import array
import sys

with open(sys.argv[1], 'rb') as in_file, open(sys.argv[2], 'w') as out_file:
    while True:
        try:
            # H = uint16_t, 16-bit unsigned integer
            a = array.array('H')
            a.fromfile(in_file, 8192)
        except EOFError:
            break
        finally:
            for x in a:
                out_file.write("{}\n".format(x))


# from collections import Counter

# c = Counter(a)
# print(c.most_common(5))

# for i in range(0, 20):
#     print(a[i])
