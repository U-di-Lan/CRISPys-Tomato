'''
The goal is to apply the crispys algorithm to gene family of solanum licopersicum (tomato).
The final output should be a library that contain sgRNA for each gene family (The sgRNA that could knockout as many members of the gene family)

First I download from PLAZA 4.0 the file of gene families (https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v4_dicots/download/index)

I used the homologus gene families.
The file name is genefamily_data.hom.csv

Lets look at him
'''

import pandas as pd


families = pd.read_csv("../genefamily_data.hom.csv", sep='\t', comment='#')

# print(families.head)

'''
The file has 3 columns:
Family name
Organism abbreviation
Gene name
'''
families.columns = ['Family', 'Organism', 'Gene']
# what are the Organisms?
# print(set(families['Organism']))

# The tomato is 'sly' so I`ll make a new file with only genes from tomato

families_tomato = families.loc[families['Organism'] == 'sly']

# print(families_tomato.head)

# how many genes?
print("Number of genes:", families_tomato.shape[0])

# #sanity
# print(set(families_tomato['Organism']))

# Write it out
families_tomato.to_csv("families_tomato.csv", sep='\t', index=False)

'''
The input to the crispys is a fasta file with exon sequence (for each gene under the same name)
So I need to extract the exon sequences.
First  I take the GFF file (also from PLAZA) and I wan to leave only the exons data
'''
pd.set_option('display.max_columns', 500)  # expand the columns display
pd.set_option('display.max_colwidth', None)  # expand the row display

gff = pd.read_csv("../annotation.all_transcripts.exon_features.sly.gff3", comment='#', sep='\t', header=None)

# print(gff.iloc[:,2]) # the elemnt columns
# Take only exons
gff_exon = gff.loc[gff.iloc[:,2] == 'exon']
#sanity
# print(set(gff_exon.iloc[:,2]))
# print(gff_exon.head)


# The gene name is in column number 8 under 'Name='
# I use regular expression to extract it
import re
# print(gff_exon.iloc[0,8])
gene_names = []
for i in range(gff_exon.shape[0]):
    name = re.findall('Name=(.+?);', gff_exon.iloc[i, 8])[0]  # take the gene name as str
    gene_names.append(name)

# print(len(gene_names))
# print(gff_exon.shape)
gff_exon['Gene_name'] = gene_names
# print(gff_exon)

# Now I have a gff of only exons
'''
The next step is to extract fasta sequences from the genome for that I`ll change the format to bed
'''
print(gff_exon)
exons_bed = gff_exon.iloc[:, [0, 3, 4, 9, 5, 6]]  # The format is: chr, start, end, name, score, strand
print(gff_exon)
# Write it out
exons_bed.to_csv("exons.bed", sep='\t', header=None, index=False)

'''
Then I use this command to create the fasta file:
 bedtools getfasta -fi sly_genome.con -bed exons.bed -s -name -fo exon_seq_fasta.fa
'''

'''
Now I have a fasta file with the exon of the genes (that belongs to the families..)
I want to edit the fasta headers so I will have only the gene name
now the header is looking like this:
 >Solyc00g005000.2::SL2.50ch00:16437-17275(+)
and I want
 >Solyc00g005000.2
'''
new_exons = open("new_exons_seq_fasta.fa", "w")
with open("exon_seq_fasta.fa", "r") as exons:
    for line in exons:
        if line.startswith(">"):
            new_exons.write((line.split(":")[0]))
            new_exons.write('\n')
        else:
            new_exons.write(line)
new_exons.close()