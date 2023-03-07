import os
import numpy as np
from Bio import SeqIO

content_liver_H3K4me1 = []
with open("../Liver/Liver-Multi-omics-data/31years/H3K4me1.bed")as f:
    for line in f:
        content_liver_H3K4me1.append(line.strip().split())    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

content_liver_H3K4me3 = []
with open("../Liver/Liver-Multi-omics-data/31years/H3K4me3.bed")as f:
    for line in f:
        content_liver_H3K4me3.append(line.strip().split())    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

content_liver_H3K9me3 = []
with open("../Liver/Liver-Multi-omics-data/31years/H3K9me3.bed")as f:
    for line in f:
        content_liver_H3K9me3.append(line.strip().split())    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

content_liver_H3K27ac = []
with open("../Liver/Liver-Multi-omics-data/31years/H3K27ac.bed")as f:
    for line in f:
        content_liver_H3K27ac.append(line.strip().split())    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

content_liver_H3K36me3 = []
with open("../Liver/Liver-Multi-omics-data/31years/H3K36me3.bed")as f:
    for line in f:
        content_liver_H3K36me3.append(line.strip().split())    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开


window_extend = [250, 400, 150, 450, 300]
featurename = '../Liver/Processed_dataset/nenhancer/'
enhancers_feature = []
H3K4me1_value_sum1 = []
H3K4me3_value_sum2 = []
H3K9me3_value_sum4 = []
H3K27ac_value_sum5 = []
H3K36me3_value_sum7 = []
for window_extend_size in window_extend:
    for record in SeqIO.parse("../Liver/Processed_dataset/nenhancer/noEnhancer_seq.fa", "fasta"):
        name = record.id
        chromosome = name.split(":")[0]
        chromosome = str(chromosome)
        # print(chromosome)
        chromosome = chromosome.strip()
        # print(chromosome1)
        enhancer_start = name.split(":")[1].split("-")[0]
        enhancer_end = name.split("-")[1]
        enhancer_start = int(enhancer_start)
        enhancer_end = int(enhancer_end)
        enhancer_center_point = int((enhancer_end - enhancer_start)/2) + enhancer_start
        start_position = enhancer_center_point - int(window_extend_size)
        end_position = enhancer_center_point + int(window_extend_size)

        if int(window_extend_size) == 250:
            value_one1 = [0]
            value_sum1 = 0
            chromStartlist1 = [0, 0]
            chromEndlist1 = [0, 0]
            chromosomelist1 = ['chr1', 'chr1']
            # content_liver_H3K4me1
            for i in content_liver_H3K4me1:
                chrom = i[0].strip()
                chrom = str(chrom)
                chromStart = int(i[1])
                chromEnd = int(i[2])
                signalValue = float(i[6])
                chromStartlist1.append(chromStart)
                chromEndlist1.append(chromEnd)
                chromosomelist1.append(chrom)
                pre_chromStart = chromStartlist1[len(chromStartlist1)-2]
                pre_chromEnd = chromEndlist1[len(chromEndlist1)-2]
                pre_chromosome = chromosomelist1[len(chromosomelist1)-2].strip()
                pre_chromStart = int(pre_chromStart)
                pre_chromEnd = int(pre_chromEnd)

                if chrom != chromosome:
                    continue
                if chromEnd <= start_position:
                    continue
                if chromStart > end_position:
                    break
                # if chrom is not pre_chromosome:
                #     continue
                if chrom is pre_chromosome and chromStart > pre_chromStart and chromEnd < pre_chromEnd:
                    continue
                if chromStart <= start_position and chromEnd > start_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one1.append(value1)
                    else:
                        base_num2 = chromEnd - start_position - 1
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one1.append(value2)
                    continue

                if chromStart > start_position and chromEnd <= end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one1.append(value1)
                    else:
                        base_num2 = chromEnd - chromStart
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one1.append(value2)
                    continue

                if chromStart <= end_position and chromEnd > end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        chromStart_new = pre_chromEnd
                        if chromStart_new > end_position:
                            break
                        else:
                            base_num = end_position - chromStart_new + 1
                            base_value = signalValue / (chromEnd - chromStart)
                            value3 = base_num * base_value
                            value_one1.append(value3)
                        continue
                    else:
                        base_num = end_position - chromStart + 1
                        base_value = signalValue / (chromEnd - chromStart)
                        value3 = base_num * base_value
                        value_one1.append(value3)
                    continue

            # content_liver_H3K4me1
            for j in value_one1:
                value_sum1 = float(j) + value_sum1
            H3K4me1_value_sum1.append(round(value_sum1, 3))

        if int(window_extend_size) == 400:
            value_one2 = [0]
            value_sum2 = 0
            chromStartlist2 = [0, 0]
            chromEndlist2 = [0, 0]
            chromosomelist2 = ['chr1', 'chr1']
            # content_liver_H3K4me3
            for i in content_liver_H3K4me3:
                chrom = i[0].strip()
                chrom = str(chrom)
                chromStart = int(i[1])
                chromEnd = int(i[2])
                signalValue = float(i[6])
                chromStartlist2.append(chromStart)
                chromEndlist2.append(chromEnd)
                chromosomelist2.append(chrom)
                pre_chromStart = chromStartlist2[len(chromStartlist2)-2]
                pre_chromEnd = chromEndlist2[len(chromEndlist2)-2]
                pre_chromosome = chromosomelist2[len(chromosomelist2)-2].strip()
                pre_chromStart = int(pre_chromStart)
                pre_chromEnd = int(pre_chromEnd)

                if chrom != chromosome:
                    continue
                if chromEnd <= start_position:
                    continue
                if chromStart > end_position:
                    break
                # if chrom is not pre_chromosome:
                #     continue
                if chrom is pre_chromosome and chromStart > pre_chromStart and chromEnd < pre_chromEnd:
                    continue
                if chromStart <= start_position and chromEnd > start_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one2.append(value1)
                    else:
                        base_num2 = chromEnd - start_position - 1
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one2.append(value2)
                    continue

                if chromStart > start_position and chromEnd <= end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one2.append(value1)
                    else:
                        base_num2 = chromEnd - chromStart
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one2.append(value2)
                    continue

                if chromStart <= end_position and chromEnd > end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        chromStart_new = pre_chromEnd
                        if chromStart_new > end_position:
                            break
                        else:
                            base_num = end_position - chromStart_new + 1
                            base_value = signalValue / (chromEnd - chromStart)
                            value3 = base_num * base_value
                            value_one2.append(value3)
                        continue
                    else:
                        base_num = end_position - chromStart + 1
                        base_value = signalValue / (chromEnd - chromStart)
                        value3 = base_num * base_value
                        value_one2.append(value3)
                    continue

            # content_liver_H3K4me3
            for j in value_one2:
                value_sum2 = float(j) + value_sum2
            H3K4me3_value_sum2.append(round(value_sum2, 3))

        if int(window_extend_size) == 150:
            start_position = start_position + 300
            end_position = end_position + 300
            value_one4 = [0]
            value_sum4 = 0
            chromStartlist4 = [0, 0]
            chromEndlist4 = [0, 0]
            chromosomelist4 = ['chr1', 'chr1']
            # content_liver_H3K9me3
            for i in content_liver_H3K9me3:
                chrom = i[0].strip()
                chrom = str(chrom)
                chromStart = int(i[1])
                chromEnd = int(i[2])
                signalValue = float(i[6])
                chromStartlist4.append(chromStart)
                chromEndlist4.append(chromEnd)
                chromosomelist4.append(chrom)
                pre_chromStart = chromStartlist4[len(chromStartlist4) - 2]
                pre_chromEnd = chromEndlist4[len(chromEndlist4) - 2]
                pre_chromosome = chromosomelist4[len(chromosomelist4) - 2].strip()
                pre_chromStart = int(pre_chromStart)
                pre_chromEnd = int(pre_chromEnd)

                if chrom != chromosome:
                    continue
                if chromEnd <= start_position:
                    continue
                if chromStart > end_position:
                    break
                # if chrom is not pre_chromosome:
                #     continue
                if chrom is pre_chromosome and chromStart > pre_chromStart and chromEnd < pre_chromEnd:
                    continue
                if chromStart <= start_position and chromEnd > start_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one4.append(value1)
                    else:
                        base_num2 = chromEnd - start_position - 1
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one4.append(value2)
                    continue

                if chromStart > start_position and chromEnd <= end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one4.append(value1)
                    else:
                        base_num2 = chromEnd - chromStart
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one4.append(value2)
                    continue

                if chromStart <= end_position and chromEnd > end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        chromStart_new = pre_chromEnd
                        if chromStart_new > end_position:
                            break
                        else:
                            base_num = end_position - chromStart_new + 1
                            base_value = signalValue / (chromEnd - chromStart)
                            value3 = base_num * base_value
                            value_one4.append(value3)
                        continue
                    else:
                        base_num = end_position - chromStart + 1
                        base_value = signalValue / (chromEnd - chromStart)
                        value3 = base_num * base_value
                        value_one4.append(value3)
                    continue

            # content_liver_H3K9me3
            for j in value_one4:
                value_sum4 = float(j) + value_sum4
            H3K9me3_value_sum4.append(round(value_sum4, 3))

        if int(window_extend_size) == 450:
            value_one5 = [0]
            value_sum5 = 0
            chromStartlist5 = [0, 0]
            chromEndlist5 = [0, 0]
            chromosomelist5 = ['chr1', 'chr1']
            # content_liver_H3K27ac
            for i in content_liver_H3K27ac:
                chrom = i[0].strip()
                chrom = str(chrom)
                chromStart = int(i[1])
                chromEnd = int(i[2])
                signalValue = float(i[6])
                chromStartlist5.append(chromStart)
                chromEndlist5.append(chromEnd)
                chromosomelist5.append(chrom)
                pre_chromStart = chromStartlist5[len(chromStartlist5) - 2]
                pre_chromEnd = chromEndlist5[len(chromEndlist5) - 2]
                pre_chromosome = chromosomelist5[len(chromosomelist5) - 2].strip()
                pre_chromStart = int(pre_chromStart)
                pre_chromEnd = int(pre_chromEnd)

                if chrom != chromosome:
                    continue
                if chromEnd <= start_position:
                    continue
                if chromStart > end_position:
                    break
                # if chrom is not pre_chromosome:
                #     continue
                if chrom is pre_chromosome and chromStart > pre_chromStart and chromEnd < pre_chromEnd:
                    continue
                if chromStart <= start_position and chromEnd > start_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one5.append(value1)
                    else:
                        base_num2 = chromEnd - start_position - 1
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one5.append(value2)
                    continue

                if chromStart > start_position and chromEnd <= end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one5.append(value1)
                    else:
                        base_num2 = chromEnd - chromStart
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one5.append(value2)
                    continue

                if chromStart <= end_position and chromEnd > end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        chromStart_new = pre_chromEnd
                        if chromStart_new > end_position:
                            break
                        else:
                            base_num = end_position - chromStart_new + 1
                            base_value = signalValue / (chromEnd - chromStart)
                            value3 = base_num * base_value
                            value_one5.append(value3)
                        continue
                    else:
                        base_num = end_position - chromStart + 1
                        base_value = signalValue / (chromEnd - chromStart)
                        value3 = base_num * base_value
                        value_one5.append(value3)
                    continue

            # content_liver_H3K27ac
            for j in value_one5:
                value_sum5 = float(j) + value_sum5
            H3K27ac_value_sum5.append(round(value_sum5, 3))

        if int(window_extend_size) == 300:
            value_one7 = [0]
            value_sum7 = 0
            chromStartlist7 = [0, 0]
            chromEndlist7 = [0, 0]
            chromosomelist7 = ['chr1', 'chr1']
            # content_liver_H3K36me3
            for i in content_liver_H3K36me3:
                chrom = i[0].strip()
                chrom = str(chrom)
                chromStart = int(i[1])
                chromEnd = int(i[2])
                signalValue = float(i[6])
                chromStartlist7.append(chromStart)
                chromEndlist7.append(chromEnd)
                chromosomelist7.append(chrom)
                pre_chromStart = chromStartlist7[len(chromStartlist7) - 2]
                pre_chromEnd = chromEndlist7[len(chromEndlist7) - 2]
                pre_chromosome = chromosomelist7[len(chromosomelist7) - 2].strip()
                pre_chromStart = int(pre_chromStart)
                pre_chromEnd = int(pre_chromEnd)

                if chrom != chromosome:
                    continue
                if chromEnd <= start_position:
                    continue
                if chromStart > end_position:
                    break
                # if chrom is not pre_chromosome:
                #     continue
                if chrom is pre_chromosome and chromStart > pre_chromStart and chromEnd < pre_chromEnd:
                    continue
                if chromStart <= start_position and chromEnd > start_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one7.append(value1)
                    else:
                        base_num2 = chromEnd - start_position - 1
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one7.append(value2)
                    continue

                if chromStart > start_position and chromEnd <= end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        base_num1 = chromEnd - pre_chromEnd
                        base_value1 = signalValue / (chromEnd - chromStart)
                        value1 = base_num1 * base_value1
                        value_one7.append(value1)
                    else:
                        base_num2 = chromEnd - chromStart
                        base_value2 = signalValue / (chromEnd - chromStart)
                        value2 = base_num2 * base_value2
                        value_one7.append(value2)
                    continue

                if chromStart <= end_position and chromEnd > end_position:
                    if chromStart <= pre_chromEnd and chrom == pre_chromosome:
                        chromStart_new = pre_chromEnd
                        if chromStart_new > end_position:
                            break
                        else:
                            base_num = end_position - chromStart_new + 1
                            base_value = signalValue / (chromEnd - chromStart)
                            value3 = base_num * base_value
                            value_one7.append(value3)
                        continue
                    else:
                        base_num = end_position - chromStart + 1
                        base_value = signalValue / (chromEnd - chromStart)
                        value3 = base_num * base_value
                        value_one7.append(value3)
                    continue

            # content_liver_H3K36me3
            for j in value_one7:
                value_sum7 = float(j) + value_sum7
            H3K36me3_value_sum7.append(round(value_sum7, 3))


        # # all_value.append(value_sum)
        # all_value = {'name': name,
        #              'content_liver_H3K4me1': value_sum1,
        #              'content_liver_H3K4me3': value_sum2,
        #              'content_liver_H3K9ac': value_sum3,
        #              'content_liver_H3K9me3': value_sum4,
        #              'content_liver_H3K27ac': value_sum5,
        #              'content_liver_H3K27me3': value_sum6,
        #              'content_liver_H3K36me3': value_sum7,
        #              'content_liver_DNase': value_sum8,
        #              'content_liver_CTCF': value_sum9}
        #
        # if value_sum1==0 and value_sum2==0 and value_sum3==0 and value_sum4==0 and value_sum5==0 and value_sum6==0 and value_sum7==0 and value_sum8==0 and value_sum9==0:
        #     all_values_is0.append(all_value)
        #
        # all_values.append(all_value)

for k in range(634):
    enhancer_feature = [H3K4me1_value_sum1[k], H3K4me3_value_sum2[k], H3K9me3_value_sum4[k], H3K27ac_value_sum5[k], H3K36me3_value_sum7[k]]
    enhancers_feature.append(enhancer_feature)

enhancers_feature = np.array(enhancers_feature)
# print(enhancers_feature)
# output the mismatch profile
np.savetxt(featurename + 'nen_feature_31years_3171.txt', enhancers_feature, encoding='utf_8_sig', delimiter=' ', fmt='%.3f')
