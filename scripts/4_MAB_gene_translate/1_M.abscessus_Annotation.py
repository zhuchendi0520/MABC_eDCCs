import sys

# 映射碱基互补
comp = {"A": "T", "T": "A", "G": "C", "C": "G"}

# 初始化数据结构
start = {}
end = {}
strand = {}
name = {}
des = {}
cat = {}
gene = []
geneha = {}
igrha = {}
code = {}
genome = ""
igr = []

# 读取基因信息文件
with open("/home/zcd/script/NTM_mapping/M.abscessus/translate/2_M.abscessus_20220725") as f1:
    lines = f1.readlines()
    for i, line in enumerate(lines):
        a = line.strip().split("\t")
        start[a[0]] = int(a[3])
        end[a[0]] = int(a[4])
        strand[a[0]] = a[2]
        name[a[0]] = a[1]
        des[a[0]] = a[5]
        cat[a[0]] = a[6]
        gene.append(a[0])
        geneha[a[0]] = 1

        # 数基因间区域
        if i > 0 and int(a[3]) > end[gene[i-1]]:
            j = f"{gene[i-1]}-{gene[i]}"
            igrha[j] = 1
            igr.append(j)
            start[j] = end[gene[i-1]] + 1
            end[j] = int(a[3]) - 1

last_igr = "EFV83_RS25425-EFV83_RS00085"
igrha[last_igr] = 1
igr.append(last_igr)
start[last_igr] = 5067231
end[last_igr] = 5068231

all_regions = gene + igr

# 读取遗传密码文件
with open("/home/zcd/script/NTM_mapping/M.abscessus/translate/3_genetic_codes") as f2:
    lines = f2.readlines()
    for line in lines:
        a = line.strip().split("\t")
        code[a[0]] = a[1]

# 读取基因组序列文件
with open("/home/zcd/script/NTM_mapping/M.abscessus.fna") as f3:
    lines = f3.readlines()
    for line in lines:
        if not line.startswith(">"):
            genome += line.strip()

k = 0  # 用来判定对该位点是否存在两个基因的重叠区域。k=0时是第一次注释，k=1是第二次注释。

# 读取突变列表文件
with open(sys.argv[1]) as f4:
    lines = f4.readlines()
    for line in lines:
        a = line.strip().split("\t")
        for i in gene:
            if start[i] <= int(a[0]) <= end[i]:  # mutation location
                if "MAB" in i:
                    print(f"{a[0]}\t{a[1]}\t{a[2]}\t-\t---\t---\t{i}\t{name[i]}\t{des[i]}\t{cat[i]}")
                elif k == 0:
                    if strand[i] == "+":
                        length = end[i] - start[i] + 1
                        seq = genome[start[i]-1:start[i]-1+length]
                        loci = int(a[0]) - start[i] + 1  # loci
                        count = loci / 3
                        ct = loci // 3
                        remain = loci % 3
                        if count == ct:
                            codon = ct
                            wd = seq[loci-3:loci]
                            mt = wd[:2] + a[2]
                        elif remain == 1:
                            codon = ct + 1
                            wd = seq[loci-1:loci+2]
                            mt = a[2] + wd[1:]
                        elif remain == 2:
                            codon = ct + 1
                            wd = seq[loci-2:loci+1]
                            mt = wd[0] + a[2] + wd[2]
                    elif strand[i] == "-":
                        length = end[i] - start[i] + 1
                        sequence = genome[start[i]-1:start[i]-1+length]
                        seq = sequence[::-1]
                        loci = end[i] - int(a[0]) + 1  # loci
                        count = loci / 3
                        ct = loci // 3
                        remain = loci % 3
                        if count == ct:
                            codon = ct
                            wd = seq[loci-3:loci]
                            mt = wd[:2] + a[2]
                        elif remain == 1:
                            codon = ct + 1
                            wd = seq[loci-1:loci+2]
                            mt = a[2] + wd[1:]
                        elif remain == 2:
                            codon = ct + 1
                            wd = seq[loci-2:loci+1]
                            mt = wd[0] + a[2] + wd[2]
                        wd = wd.translate(str.maketrans("ATGC", "TACG"))
                        mt = mt.translate(str.maketrans("ATGC", "TACG"))
                    if code.get(wd) == code.get(mt):
                        type_ = "Synonymous"
                    else:
                        type_ = "Nonsynonymous"
                    print(f"{a[0]}\t{a[1]}\t{a[2]}\t{codon}\t{type_}-{code.get(wd)}-{code.get(mt)}\t{wd}-{mt}\t{i}\t{name[i]}\t{des[i]}\t{cat[i]}")
                    type_ = ""
                if k == 1:
                    if strand[i] == "+":
                        length = end[i] - start[i] + 1
                        seq = genome[start[i]-1:start[i]-1+length]
                        loci = int(a[0]) - start[i] + 1  # loci
                        count = loci / 3
                        ct = loci // 3
                        remain = loci % 3
                        if count == ct:
                            codon = ct
                            wd = seq[loci-3:loci]
                            mt = wd[:2] + a[2]
                        elif remain == 1:
                            codon = ct + 1
                            wd = seq[loci-1:loci+2]
                            mt = a[2] + wd[1:]
                        elif remain == 2:
                            codon = ct + 1
                            wd = seq[loci-2:loci+1]
                            mt = wd[0] + a[2] + wd[2]
                    elif strand[i] == "-":
                        length = end[i] - start[i] + 1
                        sequence = genome[start[i]-1:start[i]-1+length]
                        seq = sequence[::-1]
                        loci = end[i] - int(a[0]) + 1  # loci
                        count = loci / 3
                        ct = loci // 3
                        remain = loci % 3
                        if count == ct:
                            codon = ct
                            wd = seq[loci-3:loci]
                            mt = wd[:2] + a[2]
                        elif remain == 1:
                            codon = ct + 1
                            wd = seq[loci-1:loci+2]
                            mt = a[2] + wd[1:]
                        elif remain == 2:
                            codon = ct + 1
                            wd = seq[loci-2:loci+1]
                            mt = wd[0] + a[2] + wd[2]
                        wd = wd.translate(str.maketrans("ATGC", "TACG"))
                        mt = mt.translate(str.maketrans("ATGC", "TACG"))
                    if code.get(wd) == code.get(mt):
                        type_ = "Synonymous"
                    else:
                        type_ = "Nonsynonymous"
                    print(f"{a[0]}\t{a[1]}\t{a[2]}\t{codon}\t{code.get(wd)}-{code.get(mt)}\t{wd}-{mt}\t{i}\t{name[i]}\t{des[i]}\t{cat[i]}")
                    type_ = ""
                k += 1
        k = 0
        for j in igr:
            if start[j] <= int(a[0]) <= end[j]:
                b = j.split("-")
                if j == "EFV83_RS25425-EFV83_RS00085":
                    left = int(a[0]) - start[j]
                    right = 5066832 - int(a[0])
                else:
                    left = int(a[0]) - end[b[0]]
                    right = start[b[1]] - int(a[0])
                print(f"{a[0]}\t{a[1]}\t{a[2]}\t-\t---\t{strand[b[0]]}{left}-{right}{strand[b[1]]}\t{j}\t{name[b[0]]}-{name[b[1]]}\t{des[b[0]]}##{des[b[1]]}\t{cat[b[0]]}##{cat[b[1]]}")
