"""
Auteur: Tommaso Barberis

Script permettant de modifier les positions genomiques (plus précisement les noms de contigs) d'un fichier vcf
"""
import sys

gff_file = "/data/home/tbarberis/projet_S3/localdata/ref/Drosophila-suzukii-annotation-3-ws.gff3"
gene_de = "/data/home/tbarberis/projet_S3/localdata/ref/fruit_G0.txt"

vcf_input=sys.argv[1]
vcf_output=sys.argv[2]

# parse the DE genes
with open(gene_de, "r") as de:
    lines = de.readlines()
    list_de = [x.replace("\n", "") for x in lines]

# using the gff file, crate a dict with as key the gene name and as value (chrom, start_pos)
with open(gff_file, "r")  as gff:
    lines = gff.readlines()
    lines = lines[1:]
    dict_pos = {}
    for line in lines:
        fields = line.split()
        
        if len(fields) != 1:
            if fields[2] == "gene":
                gene_name = fields[8].split(";")[0].replace("ID=", "")
                if gene_name in list_de:
                    dict_pos[gene_name] = (fields[0], fields[3]) # fields 0: chrom, fields 1: start position of the gene on the chromosome
            else:
                pass

with open(vcf_input, "rb") as vcf:
    lines = vcf.readlines()
    header = []
    new_lines = {}  
    for line in lines:
        line = str(line, encoding="latin-1") # fix parsing of the line        
        if line.startswith("#"):
            header.append(line)
        else:
            line =  line.split()
           
            gene = line[0] 
            
            if gene in list_de:
                pos_snp = line[1] # SNP position on his gene
                pos_snp = int(pos_snp) + int(dict_pos[gene][1]) - 1 # SNP position on his chromosome
                chrom = dict_pos[gene][0]
                new_line = [chrom, str(pos_snp)] + line[2:] 
                new_line = "\t".join(new_line)
                new_line = new_line + "\n"    
                if chrom in new_lines.keys(): # if the key is already in the dict or not
                    new_lines[chrom].append(new_line)
                else:
                    new_lines[chrom] = [new_line]                                
               

with open(vcf_output, "a") as new_vcf:
    for line in header:
        pass
        new_vcf.write(line)
    for chrom in new_lines.keys():        
        lines_by_chr = new_lines[chrom]

        lines_by_chr.sort(key=lambda x: int(x.split()[1])) # correctly sort snp by pos     
        for line in lines_by_chr:  
            new_vcf.write(line)
                    