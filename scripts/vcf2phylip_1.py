import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--logical", dest="log",
                  help="logical parameter", action="store_true")
parser.add_option("--param", dest="param",
                  help="numerical parameter", default=1)

(options, args) = parser.parse_args()
parser.add_option_group(group)

iupacdict = {'M': 'AC', 'R': 'AG', 'W': 'AT', 'S': 'CG', 'Y': 'CT', 'K': 'GT',
             'V': 'ACG', 'H': 'ACT', 'D': 'AGT', 'B': 'CGT', 'X': 'ACGT', 'N': 'ACGT'}

iupacrevdict = {y: x for x, y in iupacdict.items()}


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


def code(x):
    if x == "./.":
        return "N"
    elif x == "0/0":
        return REF
    elif x == "0/1":
        return iupacrevdict["".join(sorted([REF, ALT]))]
    else:
        return ALT


STRING = d(list)

X = 0
for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    if l.startswith("#"):
        a = l.rstrip().split()
        header = a[9:]
        continue

    a = l.rstrip().split()
    REF = a[3]
    ALT = a[4]
    POPS = a[9:]
    pop = [x.split(":")[0] for x in POPS]

    for i in range(len(header)):
        STRING[header[i]].append(code(pop[i]))
    X += 1

print(str(len(header)) + "\t" + str(X))
for k, v in sorted(STRING.items()):
    print(k + "\t" + "".join(v))