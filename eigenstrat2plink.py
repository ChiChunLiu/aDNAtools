#!/bin/python

from itertools import zip_longest
import numpy as np
import os
import binascii


# Eigenstrat homozygous as a0 is coded as 2
an2geno_dict = {'0': '11', '1': '10', '2': '00', '9': '01'}
sex_token = {'M': '1', 'F': '2', 'U': '0'}

def digit2string_sex(sex, token_dict):
    if sex in token_dict.keys():
        return token_dict[sex]
    else:
        return '0'

def grouper(iterable, n, fillvalue= '2'):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, '0') --> ABC DEF G00"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def geno2bed(geno, bed):
    with open(eigen_geno, 'r') as geno, open(plink_bed, 'ab') as bed:
        # plink magic header
        bed.write(b'\x6c\x1b\x01')
        
        for line in geno:
            line = grouper(line.strip(), 4)
            an_holder = []
            for an in line:
                an_holder.append(''.join(an[::-1]))

            line = ''.join(an_holder)
            line = ''.join([an2geno_dict[an] for an in list(line)])
            # Ignacio Vazquez-Abrams @ StackOverflow
            # 2072351/python-conversion-from-binary-string-to-hexadecimal
            line = '%0*X' % ((len(line) + 3) // 4, int(line, 2))
            line = binascii.unhexlify(line)

            bed.write(line)
            
def ind2fam(ind, fam):
    with open(ind, 'r') as ind, open(fam, 'a') as fam:

        if (args.population != '') & (args.ind2pop != ''):
            raise Exception('flag -p and -P cannot be used together.')
            
        nind = sum(1 for line in ind if line.strip())
        
        if args.ind2pop is not '':
            pops = []
            with open(args.ind2pop) as pfile:
                for p in pfile:
                    p = p.strip()
                    pops.append(p)
            assert len(pops) == ind
                    
        for line in ind:
            line = line.strip().split()
            line = [line[2], line[0], '0', '0', line[1], '2']
            line[4] = digit2string_sex(line[4], sex_token)

            if args.population is not '':
                line[0] = args.population
            elif args.ind2pop is not '':
                line[0] = pops.pop()

            fam.write('\t'.join(line) + '\n')


def snp2bim(snp, bim):
    with open(snp, 'r') as snp, open(bim, 'a') as bim:
        for line in snp:
            line = line.strip().split()
            line[0], line[1] = line[1], line[0]
            bim.write('\t'.join(line) + '\n')


if __name__== "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input', type = str, default = "", help = "prefix for ancestrymap input")
    parser.add_argument('-o', '--output', type = str, default = "", help = "prefix for plink output")
    parser.add_argument('-p', '--population', type = str, default = "", help = "fill in single population")
    parser.add_argument('-P', '--ind2pop', type = str, default = "", help = "list of populations")
    args = parser.parse_args()

    snp2bim(args.input + '.snp', args.output + '.bim')
    ind2fam(args.input + '.ind', args.output + '.fam')
    geno2bed(args.input + '.geno', args.output + '.bed')

