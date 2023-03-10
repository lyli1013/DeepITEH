import codecs
from keras_bert import load_trained_model_from_checkpoint, Tokenizer
import numpy as np
from keras.models import *
from keras.layers import Dense, Dropout, LSTM, GRU,Embedding, Bidirectional, Concatenate
from keras.layers import Input, Dense, Dropout, BatchNormalization, LSTM, Bidirectional
import os
import copy
import random
from keras.callbacks import ModelCheckpoint
from keras.optimizers import Adam

# Hyperparameter settings
k_merList = [2]  # kmer
kkk = k_merList[0]
EPCOHS = 10  # an arbitrary cutoff, generally defined as "one pass over the entire dataset", used to separate training into distinct phases, which is useful for logging and periodic evaluation.
BATCH_SIZE = 64  # a set of N samples. The samples in a batch are processed` independently,
# in parallel. If training, a batch results in only one update to the model.
maxlen = 718
RNN_HIDDEN_DIM = 64  # Number of neural units in the lstm layer
DROPOUT_RATIO_1 = 0.5  # proportion of neurones not used for training
DROPOUT_RATIO_2 = 0.5  # proportion of neurones not used for training
DROPOUT_RATIO_3 = 0.1
xSquare_selection = 45  # Filter out the number of SeqPose Feature

checkpoint_dir = '/data/lyli/Project/Liver_model/Seq_feature_exact/window14'  # Store model files
if not os.path.exists(checkpoint_dir):  # if the path does not exist, generate it
    os.mkdir(checkpoint_dir)

checkpoint_dir = '/data/lyli/Project/Liver_model/Seq_feature_exact/window14/checkpoints'  # Store model files
if not os.path.exists(checkpoint_dir):  # if the path does not exist, generate it
    os.mkdir(checkpoint_dir)

checkpoint_dir = '/data/lyli/Project/Liver_model/Seq_feature_exact/window14/checkpoints/final_seq_model'  # Store model files
if not os.path.exists(checkpoint_dir):  # if the path does not exist, generate it
    os.mkdir(checkpoint_dir)

output_dir = '/data/lyli/Project/Liver_model/Seq_feature_exact/window14/output'  # Store model files
if not os.path.exists(output_dir):  # if the path does not exist, generate it
    os.mkdir(output_dir)
filepath = checkpoint_dir + "/liver_seq_model_" + str(kkk) + "_mer" + ".hdf5"  # First layer model

# Read data
def read_Fasta(filename3, filename4):
    finalFasta_positive = {}
    # fileList = ['/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/25years/enhancer/train_enseq2537_300.txt']
    fileList = [filename3]
    for m in fileList:
        fasta = {}
        with open(m) as file:
            sequence = ""
            for line in file:
                if line.startswith(">"):
                    name = line[1:].rstrip()
                    fasta[name] = ''
                    continue

                fasta[name] += line.rstrip()

            finalFasta_positive = dict(finalFasta_positive, **fasta)

    finalFasta_negative = {}
    fasta = {}
    # with open('/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/25years/nenhancer/train_nenseq2537_300.txt') as file:
    with open(filename4) as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):
                name = line[1:].rstrip()
                fasta[name] = ''
                continue
            fasta[name] += line.rstrip()

        finalFasta_negative = dict(finalFasta_negative, **fasta)

    AllPositive_items = finalFasta_positive.items()
    AllNegative_items = finalFasta_negative.items()
    sequenceList = []
    for key, value in AllPositive_items:
        tmpDict = {'sequence': value, 'target': 1}  # Construct positive samples
        sequenceList.append(tmpDict)
    print("????????????", len(sequenceList))
    for key, value in AllNegative_items:
        tmpDict = {'sequence': value, 'target': 0}  # Construct negative samples
        sequenceList.append(tmpDict)
    print("????????????", len(sequenceList))
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


def Model_TransBERT(input_seq_length):
    '''
    Define the classification model
    '''
    def create_lstm(input_seq_length, rnn_hidden_dim=RNN_HIDDEN_DIM,
                    dropout_1=DROPOUT_RATIO_1, dropout_2=DROPOUT_RATIO_2, dropout_3=DROPOUT_RATIO_3):
        inputs_seq1 = Input(shape=(input_seq_length, 768))
        # inputs_omics = Input(shape=(input_omics_length,))

        bi_lstm_1 = Bidirectional(LSTM(rnn_hidden_dim, return_sequences=True))(inputs_seq1)
        bn = BatchNormalization()(bi_lstm_1)
        # attention_mul_bilstm = attention_3d_block(bn, input_length)
        drop_1 = Dropout(DROPOUT_RATIO_1)(bn)
        bi_lstm_2 = Bidirectional(LSTM(rnn_hidden_dim))(drop_1)
        drop_2 = Dropout(DROPOUT_RATIO_2)(bi_lstm_2)
        seq_dense = Dense(1, activation='sigmoid')(drop_2)

        model = Model(input=[inputs_seq1], output=[seq_dense])

        adam = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001 / EPCOHS)
        model.compile(loss='binary_crossentropy', optimizer=adam, metrics=["accuracy"])

        # dense = Dense(1, activation='sigmoid')(drop_2)
        # model = Model(inputs_seq1, inputs_omics, dense)
        # model.compile('adam', 'binary_crossentropy', metrics=['accuracy'])

        return model

    return create_lstm(input_seq_length)


# ???????????????
class OurTokenizer(Tokenizer):
  def _tokenize(self, text):
      R = []
      for i in range(len(text)-2+1):
          tmp = ''
          for j in range(i, i + 2):
              tmp += text[j]
          if tmp in self._token_dict:
              R.append(tmp)
          elif self._is_space(tmp):
              R.append('[unused100]')  # space?????????????????????[unused1]??????
          else:
              R.append('[UNK]')          # ??????????????????[UNK]
      return R

# Bert???????????????
def extract_single_features(tokenizer, text, max_length=768, max_len=512):
  # tokens = tokenizer.tokenize(text)
  indices, segments = tokenizer.encode(first=text, max_len=max_len)
  predicts = model_bert.predict([np.array([indices]), np.array([segments])])[0]

  # embeddings = []
  # for i, token in enumerate(tokens):
  #     embeddings.append(predicts[i].tolist()[:max_length])
  # embeddings = np.array(embeddings)
  return predicts


filename1 = "/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/Sequence_feature/enhancer/window14/train_enseq2537_300.txt"
filename2 = "/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/Sequence_feature/nenhancer/train_nenseq2537_300.txt"


sequenceList = read_Fasta(filename1, filename2)
random.seed(75)
random.shuffle(sequenceList)  # shuffle
# print("????????????",len(sequenceList))

wordVectorSeqList, wordVectorTarList = k_mer_split(sequenceList, k_merList)
# print(wordVectorSeqList)
# print(np.shape(wordVectorSeqList))
# print(np.shape(wordVectorSeqList[0]))
# print(type(wordVectorSeqList))
# print(type(wordVectorSeqList[0]))



config_path = '/data/lyli/Project/bert/cased_L-12_H-768_A-12/bert_config.json'
checkpoint_path = '/data/lyli/Project/bert/cased_L-12_H-768_A-12/bert_model.ckpt'
dict_path = '/data/lyli/Project/bert/cased_L-12_H-768_A-12/vocab.txt'
model_bert = load_trained_model_from_checkpoint(config_path, checkpoint_path)


# ??????Bert????????????????????????
token_dict = {}
with codecs.open(dict_path, 'r', 'utf8') as reader:
  for line in reader:
      token = line.strip()
      token_dict[token] = len(token_dict)
# print(token_dict)

# ?????????????????????
n = 1
for pairs in wordVectorSeqList:
    for token in pairs:
        if token not in token_dict:
            token_dict[token] = token_dict.pop('[unused'+str(n)+']')
            n += 1
# print(token_dict)
np.save('/data/lyli/Project/Liver_model/Seq_feature_exact/window14/ModelDict_sequence.npy', token_dict)



tokenizer = OurTokenizer(token_dict)
X_seq_train = [extract_single_features(tokenizer, x['sequence']) for x in sequenceList]
X_seq_train = np.array(X_seq_train)
Y_train = wordVectorTarList


# training Model
if __name__ == '__main__':
    model = Model_TransBERT(X_seq_train.shape[1])  # Get the first layer model
    saveBestModel = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, mode='max')
    callbacks_list = [saveBestModel]
    print('Fitting first model...')
    history = model.fit(X_seq_train, Y_train, batch_size=BATCH_SIZE, epochs=EPCOHS, callbacks=callbacks_list, validation_split=0.1, verbose=1)  # Start training


# evaluation model
filename3 = "/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/Sequence_feature/enhancer/window14/test_enseq634_300.txt"
filename4 = "/data/lyli/deep-dataset/Processed_dataset/liver_Experiment/method03/Sequence_feature/nenhancer/test_nenseq634_300.txt"
sequenceList = read_Fasta(filename3, filename4)
random.seed(75)
random.shuffle(sequenceList)  # shuffle

# wordVectorTarList is lable of test dataset
wordVectorSeqList, wordVectorTarList = k_mer_split(sequenceList, k_merList)

# load bert model
config_path = '/data/lyli/Project/bert/cased_L-12_H-768_A-12/bert_config.json'
checkpoint_path = '/data/lyli/Project/bert/cased_L-12_H-768_A-12/bert_model.ckpt'
dict_path = '/data/lyli/Project/bert/cased_L-12_H-768_A-12/vocab.txt'
model_bert = load_trained_model_from_checkpoint(config_path, checkpoint_path)

# Load SeqPoseDict
token_dict = np.load('/data/lyli/Project/Liver_model/Seq_feature_exact/window14/ModelDict_sequence.npy', allow_pickle=True).item()

# Reconstruct tokenizer
tokenizer = OurTokenizer(token_dict)

# get wordVector of test dataset
X_seq_test = [extract_single_features(tokenizer, x['sequence']) for x in sequenceList]
X_seq_test = np.array(X_seq_test)
Y_test = wordVectorTarList


# starting evaluation Model
model = Model_TransBERT(X_seq_test.shape[1])  # Get the first layer model
model.load_weights(filepath)
loss, acc = model.evaluate(X_seq_test, Y_test, batch_size=BATCH_SIZE)
print("loss:", loss)
print("acc:", acc)















