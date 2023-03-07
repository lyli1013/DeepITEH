import random
import os
import numpy as np
from Bio import SeqIO
import pysam
import math

enhancers = []
for record in SeqIO.parse("/data/lyli/Pancreatic-adenocarcinoma/Processed_dataset-new/nenhancer/02eight_cells_noenhancers1126_6001.fa", "fasta"):
    name = record.id
    enhancers.append(name)
enhancers_new = []
file = [12,22,182,206,245,307,543,583,743,835]
for i in file:
    la = int(i)
    enhancers_new.append(enhancers[la])

newfile1 = open("/data/lyli/Pancreatic-adenocarcinoma/Processed_dataset-new/nenhancer/02eight_cells_noenhancers10_6001.fa", 'w')
for i in enhancers_new:
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


