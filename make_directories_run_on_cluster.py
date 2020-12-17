'''
This script make folder for each family and in that folder write a file that contain all
the exons of genes in that family
'''



import os
import pandas as pd

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='set the max and min number of genes in family')
    parser.add_argument('--min_gene_number', '-mign', default=2, type=int)
    parser.add_argument('--max_gene_number', '-magn', default=10, type=int)
    args = parser.parse_args()

#	move_type = args.move_type
#	st = str(args.step_number)
 
 
# print(os.listdir())
# print(os.getcwd())
def make_directory_for_file (path_to_family_file, path_to_parent_directory, path_to_fasta):
    ''' The 'family_file' file has 3 columns, the first is family identifier and the third is gene name.
    This function will create a directory for each gene family and in that directory will write the fasta of all exons of all genes in that family '''

    # make a directory for all families
    if not os.path.exists(path_to_parent_directory):
        os.makedirs(path_to_parent_directory)

    # read the file with family name and gene name
    families = pd.read_csv(path_to_family_file, sep='\t')
    # take family names
    names = set(families.iloc[:, 0])

    # make a dictionary that will have the family name as key and gen names as value
    genes = {}

    for name in names:
        # create a key with the name of the family and value as empty list
        genes[name] = []

    # Take the genes belonging to that family
    for i in range(families.shape[0]):
        genes[families.iloc[i, 0]].append(families.iloc[i, 2])

    # keep only families with x number of genes
    families_limited = {}
    for family, genes_list in genes.items():
        if (len(genes[family]) <= args.max_gene_number) & (len(genes[family]) >= args.min_gene_number):
            families_limited[family] = genes_list

    for name in families_limited.keys():
        # make a folder with the name of the family
        if name not in os.listdir(path_to_parent_directory):
            os.makedirs('{}/{}'.format(path_to_parent_directory, name))

    # for each gene take its sequences
    for name in families_limited.keys():
        # go to the directory
        new_dir = path_to_parent_directory + "/" + name  # go into the new directory

        with open(path_to_fasta, "r") as fasta:  # open the fasta file with the sequences
            for line in fasta:
                if line.startswith(">"):
                    # check the name of the gene
                    gene = line[1:].strip()

                    if gene in families_limited[name]:
                        # if the gene belong to the family (i.e. in the list of genes)
                        # write it to a file
                        if os.path.isfile('{}/{}.fa'.format(new_dir, name)):  # if the file exist add th gene name and sequence
                            f = open(new_dir + "/" + name + ".txt", "a")
                            f.write(line)
                            f.write(fasta.readline())
                            f.write("\n")
                            f.close()
                        else:
                            f = open(new_dir + "/" + name + ".txt", "w")  # if the file dosent exist create it and add gene name and sequence
                            f.write(line)
                            f.write(fasta.readline())
                            f.write("\n")
                            f.close()




# test the function on fraction of the file
path_to_family_file = "/groups/itay_mayrose/udiland/crispys_tomato/families_tomato.csv"
path_to_parent_directory = "/groups/itay_mayrose/udiland/crispys_tomato/families"
path_to_fasta = r"/groups/itay_mayrose/udiland/crispys_tomato/new_exons_seq_fasta.fa"

make_directory_for_file(path_to_family_file, path_to_parent_directory, path_to_fasta)

