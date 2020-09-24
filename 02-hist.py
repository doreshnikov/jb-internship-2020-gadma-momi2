from typing import Tuple, Dict, List

import sys
from common import VCF
from pandas import DataFrame, Series, concat
import matplotlib.pyplot as plt
import re

if len(sys.argv) < 3:
    print('Please specify a path to the file and sample name')
    exit(0)
filename = sys.argv[1]
sample = sys.argv[2]

vcf = VCF.parse(filename, samples=[sample])
ploidy = vcf.ploidy_info()
assert sample in ploidy.columns and ploidy.iloc[0][sample] == 2

data = concat([
    vcf.genotype_info()[f'{sample}/GT'],
    vcf.data.loc[:, ['CHROM', 'POS']]
], axis=1)

counts = dict()
distances = []

print('Calculating distances...')
previous_chrom = ''
previous_pos = -1
for row in data.itertuples(index=False):
    if row[1] != previous_chrom:
        previous_chrom = row[1]
        print(f'done\n> chromosome {previous_chrom}: ...', end='')
        previous_pos = -1

    alleles = list(map(int, re.split(r'[|/]', row[0])))
    if alleles[0] != alleles[-1]:
        pos = row[2]
        if previous_pos != -1:
            diff = pos - previous_pos
            distances.append(diff)
            if diff not in counts.keys():
                counts[diff] = 0
            counts[diff] += 1
        previous_pos = pos

valuable = {key: value for key, value in counts.items() if value > 1}
print('done\nBuilding histogram...')
plt.hist(distances, bins=100, range=(1, 200000), density=True)
plt.show()
plt.bar(list(valuable.keys()), valuable.values())
plt.show()
