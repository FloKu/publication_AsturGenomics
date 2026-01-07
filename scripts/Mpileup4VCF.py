import sys
import gzip
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

## Viel Spass :-)

#########################################################   HELP   #########################################################################
usage="python %prog --input file --output file "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--VCF", dest="VCF", help="Input file")
parser.add_option("--Mpileup", dest="MP", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--DP", dest="DP", help="Output file",default="10")

(options, args) = parser.parse_args()
parser.add_option_group(group)

def load_data(x):
  ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
  import gzip
  if x=="-":
      y=sys.stdin
  elif x.endswith(".gz"):
      y=gzip.open(x,"rt", encoding="latin-1")
  else:
      y=open(x,"r", encoding="latin-1")
  return y

Positions=d(lambda: d(str))
out=gzip.open(options.OUT,"wt")
for l in load_data(options.VCF):
    if l.startswith("#"):
        continue
    a=l.rstrip().split()
    Chrom=a[0]
    Pos=a[1]
    Positions[Chrom][Pos]

C=1
print("reading VCF done")
for l in load_data(options.MP):
    a=l.rstrip().split()
    Chrom=a[0]
    Pos=a[1]
    if Pos in Positions[Chrom]:
        Positions[Chrom][Pos]=a[-2]
    if C%5000000==0:
        print(C,"lines Processed")
    C+=1

print("reading Mpileup done")

for l in load_data(options.VCF):
    if l.startswith("##"):
        out.write(l)
        continue
    elif l.startswith("#"):
        out.write(l.rstrip()+"\tAnisus_Mpileup\n")
        continue
    a=l.rstrip().split()
    Chrom=a[0]
    Pos=a[1]
    REF=a[3]
    ALT=a[4]

    RefCount=Positions[Chrom][Pos].count(".")+Positions[Chrom][Pos].count(",")
    AD=Positions[Chrom][Pos].upper().count(ALT)
    DP=RefCount+AD
    if DP<int(options.DP):
        GT="./."
    else:
        if RefCount>0 and AD>0:
            GT="0/1"
        elif RefCount==0 and AD>0:
            GT="1/1"
        elif RefCount>0 and AD==0:
            GT="0/0"
        else:
            GT="./."
    GQ="40"
    GL="0.0,0.0,0.0"
    if GT=="./.":
        out.write(l.rstrip()+"\t"+GT+"\n")
    else:
        out.write(l.rstrip()+"\t"+":".join([GT,str(DP),str(AD),GQ,GL])+"\n")
    C+=1
        
    
    
