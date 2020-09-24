## Common

File [common.py](common.py) contains definition of class `VCF`
designed for parsing and extracting data from .vcf files:

- `VCF.parse(filename, single_only, samples=None)` parses a given .vcf file
into pair `<metainfo, data>`
    - if `single_only` is set all non-single nucleotide polymorphisms
    are filtered out
    - if `samples` is specified, only given samples are left in data frame
- `vsf.metainfo` contains a meta-information dictionary
- `vsf.data` contains a data frame with SNP information
```
> print(vsf.data)

  CHROM      POS         ID  ...         NA00001         NA00002       NA00003
0    19      111          .  ...       0|0:10,10       0|0:10,10       0/1:3,3
1    19      112          .  ...       0|0:10,10       0|0:10,10       0/1:3,3
2    20    14370  rs6054257  ...  0|0:48:1:51,51  1|0:48:8:51,51  1/1:43:5:.,.
3    20    17330          .  ...  0|0:49:3:58,50    0|1:3:5:65,3  0/0:41:3:.,.
4    20  1110696  rs6040355  ...  1|2:21:6:23,27    2|1:2:0:18,2  2/2:35:4:.,.
5    20  1230237          .  ...  0|0:54:.:56,60  0|0:48:4:51,51  0/0:61:2:.,.
6    20  1234567  microsat1  ...         0/1:.:4        0/2:17:2      1/1:40:3
7    20  1235237          .  ...             0/0             0|0           ./.
8     X       10     rsTest  ...               0             0/1           0|2
```
- `vsf.genotype_info()` returns a data frame with information about genotypes
```
> print(vsf.genotype_info())

  NA00001/GT NA00002/GT NA00003/GT
0        0|0        0|0        0/1
1        0|0        0|0        0/1
2        0|0        1|0        1/1
3        0|0        0|1        0/0
4        1|2        2|1        2/2
5        0|0        0|0        0/0
6        0/1        0/2        1/1
7        0/0        0|0        ./.
8          0        0/1        0|2
```
- `vsf.ploidy_info()` returns a data frame with information about ploidy numbers
```
> print(vsf.ploidy_info())

        NA00001  NA00002  NA00003
ploidy        2        2        2
```

## Task 1 - ploidy

Code in [01-ploidy.py](01-ploidy.py) calculates each sample's ploidy by searching
for a maximum amount of alleles through all rows present in data frame.

## Task 2 - SNP distances histogram

Code in [02-hist.py](02-hist.py) calculates a histogram of distances between heterozygous 
alleles in one diploid organism.