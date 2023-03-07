
import numpy as np
from Bio import SeqIO

content_liver_H3K4me1 = []
with open("/data/lyli/Pancreas/pancreas-Multi-omics_data/pancreas-30years/H3K27ac-ENCFF999HBY.bed") as f:
    for line in f:
        content_liver_H3K4me1.append(line.strip().split())  # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开



featurename = '/data/lyli/Pancreas/pancreas-Multi-omics_data/pancreas-30years/Enhancer_peak_distribution/H3K27ac/'
w = 1
mult = ['window01.fa', 'window02.fa', 'window03.fa', 'window04.fa', 'window05.fa', 'window06.fa', 'window07.fa', 'window08.fa', 'window09.fa', 'window10.fa',
        'window11.fa', 'window12.fa', 'window13.fa', 'window14.fa', 'window15.fa', 'window16.fa', 'window17.fa', 'window18.fa', 'window19.fa', 'window20.fa', 'window21.fa']
# mult = ['window01.fa', 'window02.fa', 'window03.fa', 'window04.fa', 'window05.fa']
for j in mult:
    have_peak_enhancer_nums = 0
    enhancers_feature = []
    for record in SeqIO.parse("/data/lyli/Pancreas/pancreas-enhancers-windows/pancreas-enhancers-300bp/" + j, "fasta"):
        location_total = []
        distance_total = []
        signalValue_total = []
        name = record.id
        value_one1 = [0]
        chrom_enhancer = name.split(":")[0]
        chrom_enhancer = str(chrom_enhancer)
        # print(chromosome)
        chrom_enhancer = chrom_enhancer.strip()
        # print(chromosome1)
        start_enhancer = name.split(":")[1].split("-")[0]
        end_enhancer = name.split("-")[1]
        start_enhancer = int(start_enhancer)
        end_enhancer = int(end_enhancer)
        en_distance = end_enhancer - start_enhancer

        # content_liver_H3K4me1
        for i in content_liver_H3K4me1:
            chrom = i[0].strip()
            chrom = str(chrom)
            chromStart = int(i[1])
            chromEnd = int(i[2])
            signalValue = float(i[6])

            if chrom != chrom_enhancer:
                continue
            if chromEnd <= start_enhancer:
                continue
            if chromStart > end_enhancer:
                break

            if chromStart <= start_enhancer and chromEnd > start_enhancer:
                enhancer_lo = int((end_enhancer - start_enhancer) / 2) + start_enhancer
                peak_lo = int((chromEnd - chromStart) / 2) + chromStart
                location = peak_lo - enhancer_lo
                location_total.append(location)
                distancce = int(chromEnd - chromStart)
                distance_total.append(distancce)
                signalValue_total.append(signalValue)
                # have_peak_enhancer_nums += 1
                continue

            if chromStart > start_enhancer and chromEnd <= end_enhancer:
                enhancer_lo = int((end_enhancer - start_enhancer) / 2) + start_enhancer
                peak_lo = int((chromEnd - chromStart) / 2) + chromStart
                location = peak_lo - enhancer_lo
                location_total.append(location)
                distancce = int(chromEnd - chromStart)
                distance_total.append(distancce)
                signalValue_total.append(signalValue)
                # have_peak_enhancer_nums += 1
                continue

            if chromStart <= end_enhancer and chromEnd > end_enhancer:
                enhancer_lo = int((end_enhancer - start_enhancer) / 2) + start_enhancer
                peak_lo = int((chromEnd - chromStart) / 2) + chromStart
                location = peak_lo - enhancer_lo
                location_total.append(location)
                distancce = int(chromEnd - chromStart)
                distance_total.append(distancce)
                signalValue_total.append(signalValue)
                # have_peak_enhancer_nums += 1
                continue

        if len(location_total) != 0:
            have_peak_enhancer_nums += 1
        enhancer_feature = [name, en_distance, location_total, distance_total, signalValue_total]
        enhancers_feature.append(enhancer_feature)

    enhancers_feature = np.array(enhancers_feature, dtype=object)
    np.savetxt(featurename + 'window' + str(w) + '.txt', enhancers_feature, delimiter=',', fmt='%s')
    # print('非增强子组学特征不都为0的数量：', len(all_values_is0))
    # print('每种组学特征不都为0的非增强子：', all_values_is0)
    # print('.......................................................................................................................................................................................................................................................................................................')
    # print('所有的非增强子数量：', len(all_values))
    # print('所有的非增强子：', all_values)
    print('window' + str(w)+':'+str(have_peak_enhancer_nums))
    w = w + 1