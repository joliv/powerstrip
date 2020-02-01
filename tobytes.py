import array
import itertools
import sys

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)

minimum = 0
# with open(sys.argv[1], 'r') as in_file:
#     for line in in_file:
#         if int(line, 10) < minimum:
#             minimum = int(line, 10)

print("min: {}".format(minimum))

with open(sys.argv[1], 'r') as in_file, open(sys.argv[2], 'wb') as out_file:
    for chunk in grouper(in_file, 1024):
        a: array = array.array('h')
        for x in chunk:
            if x == None: continue
            a.append(int(x, 10) - minimum)
        a.tofile(out_file)
