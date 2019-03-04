#!/bin/python

import os
import numpy as np
import binascii


parser = argparse.ArgumentParser()
parser.add_argument('-b', '--input', type = str, default = "", help = "prefix for bed/bim/fam")
parser.add_argument('-o', '--output', type = str, default = "", help = "prefix for outputs")
parser.add_argument('-p', '--population', type = str, default = "", help = "one population")
parser.add_argument('-P', '--ind2pop', type = str, default = "", help = "list of populations")
args = parser.parse_args()


sex_token = {'1': 'M', '2': 'F'}
plink2eigen_geno_map = {0: 2, 1: 1, 2: np.nan, 3: 0}


def byte2binary(byte):
    byte = binascii.hexlify(byte)
    return "{:08b}".format(int(byte,16))

def digit2string_sex(sex, token_dict):
    if sex in token_dict.keys():
        return token_dict[sex]
    else:
        return 'U'

def convert_geno(bed, geno):
    with open(bed, 'rb') as bed, open(geno, 'a') as geno:
        # skip header
        bed.seek(3)
        #initialize
        counter = 1
        byte = b'\x00' # whater nonempty byte value
        binary_geno = ''

        while byte != b'':

            byte = bed.read(1)
            binary_string = byte2binary(byte)
            binary_geno += binary_string[::-1]

            if counter % nbyte_per_snp == 0:
                binary_geno = np.array(list(binary_geno)).astype(np.uint8)
                binary_geno = binary_geno[:(2 * nind)]
                gt = 2 * binary_geno[::2] + binary_geno[1::2]
                gt = np.array([plink2eigen_geno_map[g] for g in gt])
                geno.write(gt)
                binary_geno = ''

            counter += 1


def convert_snp(fam, ind):
    with open(fam, 'r') as fam, open(ind, 'a') as ind:

        if (args.population != '') & (args.ind2pop != ''):
            raise Exception('flag -p and -P cannot be used together.')

        if args.ind2pop is not '':
            pops = []
            with open(args.ind2pop) as pfile:
                for p in pfile:
                    p = p.strip()
                    pops.append(p)
        nind = 0
        for line in fam:
            nind += 1
            line = line.strip().split()
            line = [line[i] for i in [1,4,0]]
            line[1] = digit2string_sex(line[1], sex_token)

            if args.population is not '':
                line[2] = args.population
            elif args.ind2pop is not '':
                line[2] = pops.pop()

            ind.write('\t'.join(line) + '\n')
        return nind


def convert_ind(bim, snp):
    with open(bim, 'r') as bim, open(snp, 'a') as snp:
        nsnp = 0
        for line in bim:
            nsnp += 1
            line = line.strip().split()
            line[0], line[1] = line[1], line[0]
            snp.write('\t'.join(line) + '\n')

nind = convert_ind(args.input + '.bed', args.output + '.ind')
nsnp = convert_snp(args.input + '.bim', args.output + '.snp')
nbyte_per_snp = int(np.ceil(nind * 2 / 8))
convert_geno(args.input + '.fam', args.output + '.geno')

