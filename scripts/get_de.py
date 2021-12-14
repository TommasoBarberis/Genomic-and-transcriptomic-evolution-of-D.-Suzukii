"""
Ce script, à partir d'un fichier .gff et d'un fichier contenant une
liste de gènes (un gène par ligne), permet d'extraire les lignes du
.gff ayant comme valeur pour le champ 'Parent' le nom d'un des gènes
contenus dans la liste.

Auteur: Tommaso Barberis
"""
import re 

DE_gene_list = "/localdata/pandata/students/M2_projet_15/RNAseq/DE_gene/fruit_G0.txt"

with open(DE_gene_list, "r") as de_file:
    de_list = de_file.readlines()

new_list = []
for gene in de_list:
    gene = gene.replace("\n", "")
    new_list.append(gene)

de_list = new_list

gff_file = "/localdata/pandata/students/M2_projet_15/ref/exon.gff3"

new_list = []
with open(gff_file, "r") as gff:
    gff_lines = gff.readlines()
    for ind, line in enumerate(gff_lines):
        gene = line.split()[8].split(";")[1].replace("Parent=", "")
        # commnent following line to get on the 'mRNA' tag
        gene = gene.split("-")[:-2]
        gene = "-".join(gene)
        
        if gene in de_list:
            with open("DE_genes_exon_rnaseq.gff", "a") as de_gff:
                # new_list.append(gene)
                de_gff.write(line)
