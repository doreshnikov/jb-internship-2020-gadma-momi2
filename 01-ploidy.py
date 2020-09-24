import sys
from common import VCF
from pandas import DataFrame, Series
import re

if len(sys.argv) < 2:
    print('Please specify a path to the file')
    exit(0)
filename = sys.argv[1]

vcf = VCF.parse(filename, single_only=False)
ploidy = vcf.ploidy_info()

assert len(set(ploidy.loc['ploidy'])) == 1
print(ploidy.iloc[0, 0])
