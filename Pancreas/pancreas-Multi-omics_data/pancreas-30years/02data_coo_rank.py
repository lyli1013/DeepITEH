import os
import numpy as np
from Bio import SeqIO
from pandas import DataFrame

# chr1	752559	752863
# 目前这几种组学数据中的peak都不相交
omics_data = ['H3K4me1-ENCFF289OFU', 'H3K4me3-ENCFF129HJI', 'H3K9me3-ENCFF119CFU', 'H3K27ac-ENCFF999HBY', 'H3K36me3-ENCFF851NMY']
for i in omics_data:
    print("这是组学数据：", i)
    content_liver_omics = []
    with open("/data/lyli/Pancreas/pancreas-Multi-omics_data/pancreas-30years/" + i + ".bed") as f:
        for line in f:
            content_liver_omics.append(line.strip().split())  # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

    chromStartlist = [0, 0]
    chromEndlist = [0, 0]
    chromosomelist = ['chr1', 'chr1']
    for i in content_liver_omics:
        chrom = i[0].strip()
        chrom = str(chrom)
        chromStart = int(i[1])
        chromEnd = int(i[2])
        chromStartlist.append(chromStart)
        chromEndlist.append(chromEnd)
        chromosomelist.append(chrom)
        pre_chromStart = chromStartlist[len(chromStartlist) - 2]
        pre_chromEnd = chromEndlist[len(chromEndlist) - 2]
        pre_chromosome = chromosomelist[len(chromosomelist) - 2].strip()
        pre_chromStart = int(pre_chromStart)
        pre_chromEnd = int(pre_chromEnd)

        # 判断此时的起始坐标有无小于之前的起始坐标
        if chrom == pre_chromosome and chromStart < pre_chromStart:
            print("此时的起始坐标有小于之前起始坐标的peak(false)：")
            print(chrom, chromStart, chromEnd)

        # 判断当此时的和之前的起始坐标相等时，此时的终止坐标是否小于等于之前的，若是，重新排序
        if chrom == pre_chromosome and chromStart == pre_chromStart and chromEnd <= pre_chromEnd:
            print("当此时的和之前的起始坐标相等时，此时的终止坐标有小于等于之前的的peak(false)：")
            print(chrom, chromStart, chromEnd)

        # 判断当此时的起始坐标大于之前的时，此时的终止坐标是否小于或等于之前的
        if chrom == pre_chromosome and chromStart > pre_chromStart and chromEnd <= pre_chromEnd:
            print("当此时的起始坐标大于之前的时，此时的终止坐标是否小于或等于之前的peak(true)：")
            print(chrom, chromStart, chromEnd)

        # print("当此时的起始坐标大于之前的时，此时的终止坐标是否大于之前的peak(true)：")
        if chrom == pre_chromosome and chromStart > pre_chromStart and chromEnd <= pre_chromEnd:
            print("当此时的起始坐标大于之前的时，此时的终止坐标是否大于之前的peak(true)：")
            print(chrom, chromStart, chromEnd)