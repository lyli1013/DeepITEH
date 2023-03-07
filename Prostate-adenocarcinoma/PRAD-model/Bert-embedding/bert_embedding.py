import codecs
from keras_bert import load_trained_model_from_checkpoint, Tokenizer
import numpy as np
import copy
import random


# Hyperparameter settings
k_merList = [2]  # kmer
kkk = k_merList[0]


def read_Fasta(filename3, filename4):
    finalFasta_positive = {}
    fileList = [filename3]
    for m in fileList:
        fasta = {}
        with open(m) as file:
            sequence = ""
            for line in file:
                if line.startswith(">"):
                    name = line[1:].rstrip()
                    # print(name)
                    fasta[name] = ''
                    continue

                fasta[name] += line.rstrip()
            finalFasta_positive = dict(finalFasta_positive, **fasta)

    finalFasta_negative = {}
    fasta = {}
    with open(filename4) as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):
                name = line[1:].rstrip()
                fasta[name] = ''
                continue
            fasta[name] += line.rstrip()

        finalFasta_negative = dict(finalFasta_negative, **fasta)
        print(len(finalFasta_negative))
    AllPositive_items = finalFasta_positive.items()
    AllNegative_items = finalFasta_negative.items()
    sequenceList = []
    for key, value in AllPositive_items:
        tmpDict = {'sequence': value, 'target': 1}  # Construct positive samples
        sequenceList.append(tmpDict)
    print("样本个数", len(sequenceList))
    for key, value in AllNegative_items:
        tmpDict = {'sequence': value, 'target': 0}  # Construct negative samples
        sequenceList.append(tmpDict)
    print("样本个数", len(sequenceList))
    return sequenceList


# Use kmer to segment raw data
def k_mer_split(sequenceList, k_merList):
    def k_mer(sequence, k):
        # kmer algorithm
        splitSequence = []
        for i in range(len(sequence) - k + 1):
            tmp = ''
            for j in range(i, i + k):
                tmp += sequence[j]

            splitSequence.append(tmp)
        return splitSequence

    finalSequenceList = []
    for k_merNum in k_merList:
        sequenceTmpList = copy.deepcopy(sequenceList)
        for i in range(len(sequenceList)):
            sequenceTmpList[i]['sequence'] = k_mer(sequenceList[i]['sequence'], k_merNum)
            finalSequenceList.append(sequenceTmpList[i])

    # Take each word as a feature,
    avertFinal = finalSequenceList[:len(sequenceList)]
    for i in range(len(k_merList)):
        tmplist = finalSequenceList[len(sequenceList) * i:(i + 1) * len(sequenceList)]
        for dictindex in range(len(tmplist)):
            avertFinal[dictindex]['sequence'] = avertFinal[dictindex]['sequence']  # + tmplist[dictindex]['sequence']

    finalSequenceList = avertFinal  # finalSequenceList is a list, each element of which is a dictionary, and the sequence in finalSequenceList has been split by k-mer

    wordVectorSeqList = []
    wordVectorTarList = []
    for allDictSeq in finalSequenceList:
        wordVectorSeqList.append(allDictSeq['sequence'])
        wordVectorTarList.append(allDictSeq['target'])

    return wordVectorSeqList, wordVectorTarList


# 重构分词器
class OurTokenizer(Tokenizer):
  def _tokenize(self, text):
      R = []
      for i in range(len(text)-2+1):
          tmp = ''
          for j in range(i, i + 2):
              tmp += text[j]
          if tmp in self._token_dict:
              R.append(tmp)
          else:
              R.append('[UNK]')          # 剩余的字符是[UNK]
      return R


# Bert生成词向量
def extract_single_features(tokenizer, text, max_length=768, max_len=512):
  indices, segments = tokenizer.encode(first=text, max_len=max_len)
  predicts = model_bert.predict([np.array([indices]), np.array([segments])])[0]
  return predicts


filename3 = "/data/lyli/Prostate-adenocarcinoma/Processed_dataset/enhancer/train_enseq771.txt"
filename4 = "/data/lyli/Prostate-adenocarcinoma/Processed_dataset/nenhancer/train_nenseq771.txt"


sequenceList = read_Fasta(filename3, filename4)
random.seed(75)
random.shuffle(sequenceList)  # shuffle
print("样本个数",len(sequenceList))

wordVectorSeqList, wordVectorTarList = k_mer_split(sequenceList, k_merList)


config_path = '/data/lyli/Prostate-adenocarcinoma/PRAD-model/Bert-embedding/cased_L-12_H-768_A-12/bert_config.json'
checkpoint_path = '/data/lyli/Prostate-adenocarcinoma/PRAD-model/Bert-embedding/cased_L-12_H-768_A-12/bert_model.ckpt'
dict_path = '/data/lyli/Prostate-adenocarcinoma/PRAD-model/Bert-embedding/cased_L-12_H-768_A-12/vocab.txt'
model_bert = load_trained_model_from_checkpoint(config_path, checkpoint_path)


# 加载Bert词汇表，存入字典
token_dict = {}
with codecs.open(dict_path, 'r', 'utf8') as reader:
  for line in reader:
      token = line.strip()
      token_dict[token] = len(token_dict)
# print(token_dict)

# 更新词汇表字典
n = 1
for pairs in wordVectorSeqList:
    for token in pairs:
        if token not in token_dict:
            token_dict[token] = token_dict.pop('[unused'+str(n)+']')
            n += 1
# print(token_dict)
np.save('/data/lyli/Prostate-adenocarcinoma/PRAD-model/Bert-embedding/prad_seq_Dict.npy', token_dict)



tokenizer = OurTokenizer(token_dict)
X_seq_train = [extract_single_features(tokenizer, x['sequence']) for x in sequenceList]
X_seq_train = np.array(X_seq_train)
print(len(X_seq_train))
print(len(X_seq_train[0]))
print(len(X_seq_train[0][0]))
with open('/data/lyli/Prostate-adenocarcinoma/PRAD-model/Bert-embedding/PRAD_seq_train_1542_wordVector.txt', 'w') as outfile:
    for slice_2d in X_seq_train:
        np.savetxt(outfile, slice_2d, fmt='%f', delimiter=',')

