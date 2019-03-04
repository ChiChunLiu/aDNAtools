#!/bin/python
from itertools import zip_longest
import numpy as np
import binascii
import argparse
import os

# Eigenstrat homozygous as a0 is coded as 2
an2geno_dict = {'0': '11', '1': '10', '2': '00', '9': '01'}
sex_token = {'M': '1', 'F': '2', 'U': '0'}
bit_eigin2plink_dict = {'00':'11','01':'10','10':'00','11':'01'}

def digit2string_sex(sex, token_dict):
    if sex in token_dict.keys():
        return token_dict[sex]
    else:
        return '0'

# see itertools package
def grouper(iterable, n, fillvalue= '2'):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, '0') --> ABC DEF G00"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def byte2binary(byte):
    byte = binascii.hexlify(byte)
    return "{:08b}".format(int(byte,16))

def geno2bed(geno, bed):
    
    
    if geno.endswith('.geno'):
        with open(geno, 'r') as geno, open(bed, 'ab') as bed:
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
                
    elif geno.endswith('.packedancestrymapgeno'):
                
        byte = b'\x00' # whater nonempty byte value
        nbyte_per_snp = max(48,int(np.ceil(nind*2/8)))
        padding = 8 - (nind*2) % 8
        
        with open(geno, 'rb') as geno, open(bed, 'ab') as bed:
            #skip header
            geno.seek(nbyte_per_snp)
            # write plink magic header
            bed.write(b'\x6c\x1b\x01')
            
            counter = 1
            while byte != b'':

                geno_byte = geno.read(1)
                if geno_byte == b'':
                    break
                binary_string = byte2binary(geno_byte)
                binary_geno = list(grouper(binary_string, 2))
                # plink and ancestrymap fill bytes right-to-left and left-to-right
                binary_geno = [bit_eigin2plink_dict[''.join(g)] for g in binary_geno[::-1]]
                binary_plink = ''.join(binary_geno)
                # map paddings if extra bits
                if (counter % nbyte_per_snp == 0) & (padding > 0):
                    binary_plink = '0' * padding + binary_plink[padding:]
                counter += 1
                
                line = '%0*X' % ((len(binary_plink) + 3) // 4, int(binary_plink, 2))
                line = binascii.unhexlify(line)
                bed.write(line)
    else:
        raise Exception('input should be in either packed or unpacked ancestry map format.')
            
def ind2fam(ind, fam):
    with open(ind, 'r') as ind, open(fam, 'a') as fam:

        if (args.population != '') & (args.ind2pop != ''):
            raise Exception('flag -p and -P cannot be used together.')
        
        if args.ind2pop is not '':
            pops = []
            with open(args.ind2pop) as pfile:
                for p in pfile:
                    p = p.strip()
                    pops.append(p)
            assert len(pops) == ind
        
        nind = 0
        for line in ind:
            nind += 1
            line = line.strip().split()
            line = [line[2], line[0], '0', '0', line[1], '2']
            line[4] = digit2string_sex(line[4], sex_token)

            if args.population is not '':
                line[0] = args.population
            elif args.ind2pop is not '':
                line[0] = pops.pop()  
            fam.write('\t'.join(line) + '\n')
        return nind

def snp2bim(snp, bim):
    with open(snp, 'r') as snp, open(bim, 'a') as bim:
        for line in snp:
            line = line.strip().split()
            line[0], line[1] = line[1], line[0]
            bim.write('\t'.join(line) + '\n')


if __name__== "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type = str, default = "", help = "prefix for ancestrymap input")
    parser.add_argument('-c', '--compressed', action='store_true', help = "input geno packed or not")
    parser.add_argument('-o', '--output', type = str, default = "", help = "prefix for plink output")
    parser.add_argument('-p', '--population', type = str, default = "", help = "fill in single population")
    parser.add_argument('-P', '--ind2pop', type = str, default = "", help = "list of populations")
    args = parser.parse_args()
    
    snp2bim(args.input + '.snp', args.output + '.bim')
    nind = ind2fam(args.input + '.ind', args.output + '.fam')
    if args.compressed == False:
        geno2bed(args.input + '.geno', args.output + '.bed')
    else:
        geno2bed(args.input + '.packedancestrymapgeno', args.output + '.bed')
