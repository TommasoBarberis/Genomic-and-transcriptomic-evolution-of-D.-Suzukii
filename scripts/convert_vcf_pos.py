import sys

gff_file = "/data/home/tbarberis/projet_S3/localdata/ref/Drosophila-suzukii-annotation-3-ws.gff3"
gene_de = "/data/home/tbarberis/projet_S3/localdata/ref/fruit_G0.txt"

vcf_input=sys.argv[1]
vcf_output=sys.argv[2]

with open(gene_de, "r") as de:
    lines = de.readlines()
    list_de = [x.replace("\n", "") for x in lines]

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
                    dict_pos[gene_name] = (fields[0], fields[3])
            else:
                pass

with open(vcf_input, "rb") as vcf:
    lines = vcf.readlines()
    header = []
    new_lines = []    
    for line in lines:
        line = str(line, encoding="latin-1") # fix parsing of the line        
        if line.startswith("#"):
            header.append(line)
        else:
            line =  line.split()
            # line = [str(x, encoding="utf-8") for x in line]
            gene = line[0] 
            
            if gene in list_de:
                start = line[1] 
                start = int(start) + int(dict_pos[gene][1])
                new_line = [dict_pos[gene][0], str(start)] + line[2:] 
                new_line = "\t".join(new_line)
                new_line = new_line + "\n"                                    
                new_lines.append(new_line)


with open(vcf_output, "a") as new_vcf:
    for line in header:
        new_vcf.write(line)
    for new_line in new_lines:
        new_vcf.write(new_line)
                    