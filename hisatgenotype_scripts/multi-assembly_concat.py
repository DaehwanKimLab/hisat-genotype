#!/usr/bin/env python3

import sys, os, subprocess, re
import glob

def concat_f(fns):
    data = {}
    genes = []
    for fn in fns:
        name = fn.split("-hla-")[1]
        data[name] = {}
        with open(fn, "r") as ifi:
            gene, group = '', []
            for line in ifi:
                line = line.strip()
                if "hisat2" in line or not line:
                    continue
                if "aligned" in line:
                    if gene and group:
                        if gene not in genes:
                            genes.append(gene)

                        data[name][gene] = group

                    gene, group = '', []
                    group.append(line)
                    continue

                if not gene:
                    assert line.startswith("1")
                    gene = line.split()[1].split('*')[0]

                group.append(line)

            if gene and group:
                if gene not in genes:
                    genes.append(gene)
                data[name][gene] = group

    return data, genes

def format_f(data,genes):
    matrix = [[] for x in range((len(genes)*21+1))]
    sort_genes = sorted(genes)

    for name, gene_dic in data.items():
        matrix[0].append(name)

        ptr = 1
        for gene in sort_genes:
            if gene not in gene_dic:
                reform_group = ["" for x in range(11)]
            else:
                reform_group = []
                ranked = False
                for x in data[name][gene]:
                    if "aligned" in x:
                        reform_group.append(x)
                        continue

                    if not ranked and "ranked" in x:
                        if len(reform_group) < 11:
                            for i in range(11 - len(reform_group)):
                                reform_group.append('')

                    reform_group.append(x)
            if len(reform_group) < 21:
                for i in range(21-len(reform_group)):
                    reform_group.append('')

            for i in range(len(reform_group)):
                matrix[i+ptr].append(reform_group[i])
            ptr += 21

    return matrix

def write_matrix(matrix):
    for i in range(len(matrix)):
        line = ",".join(matrix[i])
        print(line)

if __name__ == '__main__':
    fns = glob.glob('*.report')
    data, genes = concat_f(fns)
    matrix = format_f(data,genes)
    write_matrix(matrix)


