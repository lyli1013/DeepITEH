import numpy as np

featurename = '/data/lyli/Lung/Lung-Multi-omics_data/lower_lobe_of_left_lung_59years_hg38/4line-hg38/'

lung_omics_datas = ['H3K4me1-ENCFF479WUX', 'H3K4me3-ENCFF180MOQ', 'H3K9me3-ENCFF180RTK', 'H3K27ac-ENCFF011IHJ', 'H3K36me3-ENCFF769ZRI']
for i in lung_omics_datas:
    lung_omics_data = []
    with open("/data/lyli/Lung/Lung-Multi-omics_data/lower_lobe_of_left_lung_59years_hg38/"+i+".bed")as f:
        for line in f:
            lung_omics_data.append(line.strip().split())    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

    peaks = []
    for j in lung_omics_data:
        chrom = j[0].strip()
        chrom = str(chrom)
        chromStart = int(j[1])
        chromEnd = int(j[2])
        signalValue = float(j[6])
        peak = [chrom, chromStart, chromEnd, signalValue]
        peaks.append(peak)

    peaks = np.array(peaks)
    np.savetxt(featurename + i + '.bed', peaks, delimiter=' ', fmt='%s')