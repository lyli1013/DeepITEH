import random
import numpy as np
from Bio import SeqIO
import pysam


enhancers = []
for record in SeqIO.parse("/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/31years/enhancer/enhancer_seq_3171_300.fa", "fasta"):
    name = record.id
    enhancers.append(name)

en_feature = []
with open("/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/31years/enhancer/Enhancer_OmicsFeature_RE_AE.txt", encoding='utf_8_sig')as f:
    for line in f:
        en_feature.append(line.strip().split())    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

newfile1 = open("/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/31years/enhancer/train_enseq2537_300.txt", 'w')
random.seed(1)
test_enhancers = random.sample(enhancers, 634)
for i in test_enhancers:
    for j in enhancers:
        if i == j:
            enhancers.remove(j)
            break
for i in enhancers:
    chromosome = i.split(":")[0]
    chromosome = str(chromosome)
    # print(chromosome)
    chromosome = chromosome.strip()
    # print(chromosome1)
    start_position = i.split(":")[1].split("-")[0]
    end_position = i.split("-")[1]
    start_position = int(start_position)
    end_position = int(end_position)
    ref = '/data/lyli/Hg19/tolhg19.fa'
    with pysam.FastaFile(ref) as fa_obj:
        seq = fa_obj.fetch(chromosome, start_position, end_position)

    newfile1.write(">" + chromosome + ":" + str(start_position) + "-" + str(end_position) + "\n" + seq + "\n")
newfile1.close()

newfile2 = open("/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/31years/enhancer/test_enseq634_300.txt", 'w')
for i in test_enhancers:
    chromosome = i.split(":")[0]
    chromosome = str(chromosome)
    # print(chromosome)
    chromosome = chromosome.strip()
    # print(chromosome1)
    start_position = i.split(":")[1].split("-")[0]
    end_position = i.split("-")[1]
    start_position = int(start_position)
    end_position = int(end_position)
    ref = '/data/lyli/Hg19/tolhg19.fa'
    with pysam.FastaFile(ref) as fa_obj:
        seq = fa_obj.fetch(chromosome, start_position, end_position)

    newfile2.write(">" + chromosome + ":" + str(start_position) + "-" + str(end_position) + "\n" + seq + "\n")
newfile2.close()

random.seed(1)
test_feature = random.sample(en_feature, 634)
for i in test_feature:
    for j in en_feature:
        if i == j:
            en_feature.remove(j)
            break

test_feature = np.array(test_feature)
en_feature = np.array(en_feature)
np.savetxt("/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/31years/enhancer/test_enfeature634.txt", test_feature, encoding='utf_8_sig', fmt='%s', delimiter=' ')
np.savetxt("/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/31years/enhancer/train_enfeature2537.txt", en_feature, encoding='utf_8_sig', fmt='%s', delimiter=' ')

