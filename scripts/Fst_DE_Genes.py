import statistics

def FstGenesDE(file1,file2):
    with open(file1, "r") as f:
        lst=[]
        DE_lst=[]
        for l in f:
            lst = l.split(" ")
            lst[2]= lst[2].replace("\n", "")
            lst[1] = int(lst[1])
            lst[2] = int(lst[2])
            DE_lst.append(lst)

    with open(file2, "r") as f:
        n=0
        lstFst=[]
        for l in f:
            n+=1
            lst1 = l.split(" ")
            lst1[1]=int(lst1[1])
            lst1[2]= lst1[2].replace("\n", "")
            lstFst.append(lst1)

    chro =[]
    for i in DE_lst:
        chro.append(i[0])
    uniqChro = set(chro)

    dicoFst= {}
    lstTemp =[]
    for l in lstFst:
        for ch in uniqChro:
            if l[0]==ch:
                lstTemp.append(l)
                dicoFst[ch]=lstTemp

    def dichotomie(t, v):
        a = 0
        b = len(t) - 1
        while a <= b:
            m = (a + b) // 2
            if int(t[m][1]) == v:
                return t[m]
            elif int(t[m][1]) < v:
                a = m + 1
            else:
                b = m - 1

    geneFst = []
    for ch in uniqChro:
        for l in DE_lst:
            if l[0] == ch:
                for g in range(l[1],l[2]+1):
                    if ch in dicoFst:
                        if dichotomie(dicoFst[ch],g):
                            geneFst.append(dichotomie(dicoFst[ch],g))

    som=[]
    for i in geneFst:
        som.append(float(i[2]))
    print(f'Nombres de SNPs:{len(som)} \nFst moyen: {statistics.mean(som)} \nFst écart type: {statistics.pstdev(som)}')
    return geneFst
PosDE = "/home/bertrand/Cours Bioinfo/Projet_3/pos_DE.txt"
FstFile = "/home/bertrand/Cours Bioinfo/Projet_3/Fst_1-2.txt"

FstGenesDE(PosDE,FstFile)

