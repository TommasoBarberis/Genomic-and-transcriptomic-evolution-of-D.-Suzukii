import sys

gff_file = "/data/home/tbarberis/projet_S3/localdata/ref/Drosophila-suzukii-annotation-3-ws.gff3"

input_vcf = sys.argv[1]

with open(gff_file, "r")  as gff:
    lines = gff.readlines()
    lines = lines[1:]
    dict_gene = {}
    for line in lines:
        fields = line.split()
        
        if len(fields) != 1:
            if fields[2] == "gene":
                gene_name = fields[8].split(";")[0].replace("ID=", "")
                chrom = fields[0]
                start = fields[3]
                end = fields[4]
                if chrom in dict_gene.keys():
                    dict_gene[chrom].append((gene_name, start, end))
                else:
                    dict_gene[chrom] = [(gene_name, start, end)]
            else:
                pass

c = 0 # counter

with open(input_vcf, "rb") as vcf_file:
    lines = vcf_file.readlines()
    for line in lines:
        if line.startswith(b"#"):
            pass 
        else:
            line =  line.split()
            line = [str(x, encoding="utf-8") for x in line]
            chrom = line[0]
            pos = line[1]

            try:
                gene_list = dict_gene[chrom]
                new_list = gene_list
                for gene in gene_list:
                    if pos >= gene[1] and pos <= gene[2]:
                        c += 1
                        new_list.remove(gene)
                        dict_gene[chrom] = new_list
                        break
                        
            except:
                pass
            

print(c)            