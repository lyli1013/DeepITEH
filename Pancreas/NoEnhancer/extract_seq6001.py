from Bio import SeqIO
import pysam

no_enhancers = []
for record in SeqIO.parse("/data/lyli/Pancreas/NoEnhancer/01eight_cells_noenhancers3171.fa", "fasta"):
    name = record.id
    no_enhancers.append(name)

newfile = open("/data/lyli/Pancreas/NoEnhancer/01eight_cells_noenhancers3171_6001.fa", 'w')


for i in no_enhancers:
    chromosome = i.split(":")[0]
    chromosome = str(chromosome)
    # print(chromosome)
    chromosome = chromosome.strip()
    # print(chromosome1)
    start_position = i.split(":")[1].split("-")[0]
    end_position = i.split("-")[1]
    start_position = int(start_position)
    end_position = int(end_position)
    enhancer_center_point = int((end_position - start_position) / 2) + 1 + start_position
    start_position = enhancer_center_point - 3000
    end_position = enhancer_center_point + 3001
    ref = '/data/lyli/Hg19/tolhg19.fa'
    with pysam.FastaFile(ref) as fa_obj:
        seq = fa_obj.fetch(chromosome, start_position, end_position)

    newfile.write(">" + chromosome + ":" + str(start_position) + "-" + str(end_position) + "\n" + seq + "\n")
newfile.close()