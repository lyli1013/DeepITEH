import random
from Bio import SeqIO
import pysam

windows = ['window01', 'window02', 'window03', 'window04', 'window05', 'window06', 'window07', 'window08', 'window09',
           'window10',
           'window11', 'window12', 'window13', 'window14', 'window15', 'window16', 'window17', 'window18', 'window19',
           'window20', 'window21']

for window in windows:
    enhancers = []
    for record in SeqIO.parse("/data/lyli/Liver-hepatocellular-carcinoma/LIHC-enhancers-windows/LIHC-enhancers-300bp/" + window + ".fa", "fasta"):
        name = record.id
        enhancers.append(name)

    newfile1 = open("/data/lyli/Liver-hepatocellular-carcinoma/LIHC-enhancers-windows/LIHC-enhancers-300bp/" + window + "_train_354.txt", 'w')
    random.seed(1)
    test_enhancers = random.sample(enhancers, 88)
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

    newfile2 = open("/data/lyli/Liver-hepatocellular-carcinoma/LIHC-enhancers-windows/LIHC-enhancers-300bp/" + window + "_test_88.txt", 'w')
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