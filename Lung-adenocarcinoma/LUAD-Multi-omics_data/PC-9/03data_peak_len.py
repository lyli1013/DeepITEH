import os
import numpy as np
from Bio import SeqIO
from pandas import DataFrame

# chr1	752559	752863

omics_data = ['H3K4me1-ENCFF487TPY', 'H3K4me3-ENCFF547MFY', 'H3K9me3-ENCFF657HCS', 'H3K27ac-ENCFF735GKU', 'H3K36me3-ENCFF896YKX']
for i in omics_data:
    print("这是组学数据：", i)
    content_liver_omics = []
    with open("/data/lyli/Lung-adenocarcinoma/LUAD-Multi-omics_data/PC-9/" + i + ".bed")as f:
        for line in f:
            content_liver_omics.append(line.strip().split())    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

    n = 0
    sum = 0
    omics_len = []
    for j in content_liver_omics:
        chrom = j[0].strip()
        chrom = str(chrom)
        chromStart = int(j[1])
        chromEnd = int(j[2])
        len = chromEnd - chromStart
        omics_len.append(len)
        sum = sum + len
        n = n + 1
    print("，有"+str(n)+"个peak。\n"+"peak的总长度为："+str(sum) +"\n"+ "peak的平均长度为："+str(sum/n))
    newfile = open("/data/lyli/Lung-adenocarcinoma/LUAD-Multi-omics_data/PC-9/omics_data_len/" + i + ".txt", 'w')
    for a in omics_len:
        newfile.write(str(a) + "\n")
    newfile.close()