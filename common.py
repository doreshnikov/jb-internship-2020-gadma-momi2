from typing import List, Dict, Union, Optional, Tuple, Type

from pandas import DataFrame, Series, concat
import re


class VCF:
    keywords = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    meta_pattern = re.compile(r'(?P<key>\w+)=(?P<value>"[\w ]+"|\w+)')

    MetainfoType = Dict[str, List[Union[str, Dict[str, str]]]]

    def __init__(self, metainfo: MetainfoType, data: DataFrame):
        self.metainfo = metainfo
        self.data = data

        self.__genotype = None
        self.__ploidy = None

    @staticmethod
    def parse(filename: str, single_only: bool = True, samples: Union[List[str], None] = None) -> 'VCF':
        print(f'Parsing {filename}...')

        lines = list(map(str.strip, open(filename).readlines()))
        metainfo: MetainfoType = dict()
        header: List[str] = []

        counter, block = 0, len(lines) // 10
        tmp_data = []
        for line in lines:
            counter += 1
            if line.startswith('##'):
                key, value = line[2:].split('=', 1)
                if key not in metainfo.keys():
                    metainfo[key] = []
                if value.startswith('<'):
                    metainfo[key].append(dict())
                    for group in re.finditer(VCF.meta_pattern, value):
                        metainfo[key][-1][group['key']] = group['value'].strip('"')
                else:
                    metainfo[key].append(value)
            elif line.startswith('#'):
                print('> metainfo read')
                header = line[1:].split()
                print('> header read\n> ... %: ', end='')
            else:
                if counter % block == 0:
                    print(f'{counter // block * 10}%', end=' ')
                row: List[Any] = line.split()
                assert len(row) == len(header)
                row[1] = int(row[1])
                tmp_data.append(row)

        print('\n> data read')
        data = DataFrame(tmp_data, columns=header)
        print('> data frame created')
        if single_only:
            data = data[
                data['REF'].apply(len) == 1 |
                data['ALT'].apply(lambda alt: all(map(lambda n: len(n) == 1, alt.split(','))))
                ]
            print('> non-single nps removed')

        vcf = VCF(metainfo, data.reset_index(drop=True))
        if samples is not None:
            print('Removing unneeded samples...')
            all_samples = vcf.sample_names()
            vcf.data = vcf.data.drop(columns=list(filter(lambda c: c not in samples, all_samples)))
            print('> done')

        return vcf

    def sample_names(self) -> List[str]:
        return list(filter(lambda c: c not in VCF.keywords, self.data.columns))

    def genotype_info(self) -> DataFrame:
        if self.__genotype is not None:
            return self.__genotype

        print('Extracting genotypes...')
        info = DataFrame()

        genotype_id = 'GT'
        candidates = list(filter(
            lambda f: isinstance(f, dict) and f.get('Description', '') == 'Genotype',
            self.metainfo.get('FORMAT', [])
        ))
        if len(candidates) > 0:
            assert len(candidates) == 1
            genotype_id = candidates[0]['ID']
        print(f'> genotype FORMAT ID: {genotype_id}')

        info['$idx'] = self.data['FORMAT'].apply(lambda f: f.split(':').index(genotype_id))
        print('> genotype index in format calculated')
        for sample in self.sample_names():
            # noinspection PyTypeChecker
            info[f'{sample}/GT'] = concat([self.data[sample], info['$idx']], axis=1).apply(
                lambda f: f[0].split(':')[f[1]],
                axis=1
            )
            print(f'> sample {sample} done')

        self.__genotype = info.drop(columns=['$idx']).reset_index(drop=True)
        print('> genotypes extracted')
        return self.__genotype

    def ploidy_info(self) -> DataFrame:
        if self.__ploidy is not None:
            return self.__ploidy

        genotype_info = self.genotype_info()
        print('Calculating ploidy based on genotypes...')
        info = DataFrame(index=['ploidy'])

        for column in genotype_info.columns:
            sample = column[:-3]
            info[sample] = [Series.max(genotype_info[column].apply(
                lambda gt: len(re.split(r'[|/]', gt))
            ))]
            print(f'> sample {sample} done')

        self.__ploidy = info
        print('> ploidy calculated')
        return self.__ploidy


if __name__ == '__main__':
    example = VCF.parse('resources/example-01.vcf', single_only=False)
    print(example.data)
    print(example.genotype_info())
    print(example.ploidy_info())
