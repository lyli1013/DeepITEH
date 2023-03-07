from Bio import SeqIO
import pysam

enhancers = []
for record in SeqIO.parse("/data/lyli/Pancreas/NoEnhancer/NoEnhancer3368_6001.fa", "fasta"):
    name = record.id
    enhancers.append(name)

newfile01 = open("/data/lyli/Pancreas/NoEnhancer/window11.fa", 'w')

# 通过原始序列6001bp得到中心碱基，从中心碱基向上下游各扩展150bp，得到增强子序列中间300bp大小的窗口
for enhancer in enhancers:
    chromosome = enhancer.split(":")[0]
    chromosome = str(chromosome)
    # print(chromosome)
    chromosome = chromosome.strip()
    # print(chromosome1)
    start_position = enhancer.split(":")[1].split("-")[0]
    end_position = enhancer.split("-")[1]
    start_position = int(start_position)
    end_position = int(end_position)
    enhancer_center_point = int((end_position - start_position) / 2) + start_position
    start_position = enhancer_center_point - 150
    end_position = enhancer_center_point + 150
    ref = '/data/lyli/Hg19/tolhg19.fa'
    with pysam.FastaFile(ref) as fa_obj:
        seq = fa_obj.fetch(chromosome, start_position, end_position)

    newfile01.write(">" + chromosome + ":" + str(start_position) + "-" + str(end_position) + "\n" + seq + "\n")
newfile01.close()

# 从中间300bp的窗口往上游继续取大小为300bp的窗口，直至把增强子全部划分
left = ['window10.fa', 'window09.fa', 'window08.fa', 'window07.fa', 'window06.fa', 'window05.fa', 'window04.fa', 'window03.fa', 'window02.fa', 'window01.fa']
for i in ['window11.fa', 'window10.fa', 'window09.fa', 'window08.fa', 'window07.fa', 'window06.fa', 'window05.fa', 'window04.fa', 'window03.fa', 'window02.fa']:
    enhancers = []
    for record in SeqIO.parse("/data/lyli/Pancreas/NoEnhancer/" + i, "fasta"):
        name = record.id
        enhancers.append(name)
    for j in left:
        newfile02 = open("/data/lyli/Pancreas/NoEnhancer/" + j, 'w')
        for enhancer in enhancers:
            chromosome = enhancer.split(":")[0]
            chromosome = str(chromosome)
            # print(chromosome)
            chromosome = chromosome.strip()
            # print(chromosome1)
            start_position = enhancer.split(":")[1].split("-")[0]
            end_position = int(start_position)
            start_position = int(end_position) - 300
            ref = '/data/lyli/Hg19/tolhg19.fa'
            with pysam.FastaFile(ref) as fa_obj:
                seq = fa_obj.fetch(chromosome, start_position, end_position)

            newfile02.write(">" + chromosome + ":" + str(start_position) + "-" + str(end_position) + "\n" + seq + "\n")
        newfile02.close()
        left.remove(j)
        break

# 从中间300bp的窗口往下游继续取大小为300bp的窗口，直至把增强子全部划分
right = ['window12.fa', 'window13.fa', 'window14.fa', 'window15.fa', 'window16.fa', 'window17.fa', 'window18.fa', 'window19.fa', 'window20.fa', 'window21.fa']
for i in ['window11.fa', 'window12.fa', 'window13.fa', 'window14.fa', 'window15.fa', 'window16.fa', 'window17.fa', 'window18.fa', 'window19.fa', 'window20.fa']:
    enhancers = []
    for record in SeqIO.parse("/data/lyli/Pancreas/NoEnhancer/" + i, "fasta"):
        name = record.id
        enhancers.append(name)
    for j in right:
        newfile03 = open("/data/lyli/Pancreas/NoEnhancer/" + j, 'w')
        for enhancer in enhancers:
            chromosome = enhancer.split(":")[0]
            chromosome = str(chromosome)
            # print(chromosome)
            chromosome = chromosome.strip()
            # print(chromosome1)
            end_position = enhancer.split("-")[1]
            start_position = int(end_position)
            end_position = int(start_position) + 300
            ref = '/data/lyli/Hg19/tolhg19.fa'
            with pysam.FastaFile(ref) as fa_obj:
                seq = fa_obj.fetch(chromosome, start_position, end_position)

            newfile03.write(">" + chromosome + ":" + str(start_position) + "-" + str(end_position) + "\n" + seq + "\n")
        newfile03.close()
        right.remove(j)
        break
