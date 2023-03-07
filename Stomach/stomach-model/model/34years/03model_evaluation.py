import pandas as pd
from keras_bert import load_trained_model_from_checkpoint, Tokenizer
import numpy as np
from keras.models import *
from keras.layers import Dense, Dropout, LSTM, GRU, Embedding, Bidirectional, Concatenate
from keras.layers import Input, Dense, Dropout, BatchNormalization, LSTM, Bidirectional
import os
import copy
import random
from keras.optimizers import Adam
from sklearn.metrics import roc_auc_score

# Hyperparameter settings
k_merList = [2]  # kmer
kkk = k_merList[0]
BATCH_SIZE = 64  # best_parameter
EPCOHS = 10    # best_parameter
RNN_HIDDEN_DIM = 64  # Number of neural units in the lstm layer
DROPOUT_RATIO_1 = 0.5  # proportion of neurones not used for training
DROPOUT_RATIO_2 = 0.5  # proportion of neurones not used for training
DROPOUT_RATIO_3 = 0.1


checkpoint_dir = '/data/lyli/Stomach/stomach-model/Bert_LSTM_DNN_model/34years/final_model_34years'  # Store model files
if not os.path.exists(checkpoint_dir):  # if the path does not exist, generate it
    os.mkdir(checkpoint_dir)

model_weight_filepath = checkpoint_dir + "/stomach_34years_model.h5"  # First layer model
model_structure_filepath = checkpoint_dir + "/stomach_34years_model.json"  # First layer model



# Read data
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

# multi_omics data features
def load_data(filename1, filename2):
    filelist = [filename1, filename2]
    dataMat = []
    for file in filelist:
        file0 = np.loadtxt(file, encoding='utf_8_sig', delimiter=' ')
        for line in file0:
            curLine = list(line)
            floatLine = map(float, curLine)  # 这里使用的是map函数直接把数据转化成为float类型
            floatLine = list(floatLine)
            dataMat.append(floatLine[0:7])


    random.seed(75)
    random.shuffle(dataMat)  # shuffle

    labelMat = []
    for i in dataMat:
        labelMat.append(int(i[-1]))
        del (i[-1])

    dataMat = np.array(dataMat)
    labelMat = np.array(labelMat)
    return dataMat, labelMat


def Model_TransBERT(input_seq_length, input_omics_length):
    '''
    Define the classification model
    '''
    def create_lstm(input_seq_length, input_omics_length, rnn_hidden_dim=RNN_HIDDEN_DIM,
                    dropout_1=DROPOUT_RATIO_1, dropout_2=DROPOUT_RATIO_2, dropout_3=DROPOUT_RATIO_3):
        inputs_seq1 = Input(shape=(input_seq_length, 768))
        inputs_omics = Input(shape=(input_omics_length,))

        bi_lstm_1 = Bidirectional(LSTM(rnn_hidden_dim, return_sequences=True))(inputs_seq1)
        bn = BatchNormalization()(bi_lstm_1)
        drop_1 = Dropout(DROPOUT_RATIO_1)(bn)
        bi_lstm_2 = Bidirectional(LSTM(rnn_hidden_dim))(drop_1)
        drop_2 = Dropout(DROPOUT_RATIO_2)(bi_lstm_2)
        seq_dense = Dense(1, activation='sigmoid')(drop_2)
        print("drop_2格式")
        print(np.shape(drop_2))

        omics_dense1 = Dense(50, name="chip_first_layer")(inputs_omics)
        omics_dense2 = Dropout(DROPOUT_RATIO_3)(omics_dense1)
        omics_dense3 = Dense(4, name="chip_second_layer", activation="relu")(omics_dense2)
        omics_dense4 = BatchNormalization()(omics_dense3)


        print(np.shape(omics_dense4))
        print("omics_dense4格式")
        x = Concatenate(axis=1)([seq_dense, omics_dense4])

        pred_output = Dense(1, activation='sigmoid')(x)
        model = Model(input=[inputs_seq1, inputs_omics], output=[pred_output])

        adam = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001 / EPCOHS)
        model.compile(loss='binary_crossentropy', optimizer=adam, metrics=["accuracy"])

        return model

    return create_lstm(input_seq_length, input_omics_length)


def format_predict(y_pred_prob):
    y_pre = []
    for prob_ele in y_pred_prob:
        if prob_ele[0] >= 0.5:  # class=1
            y_pre.append(1)
        else:
            y_pre.append(0)

    return y_pre


# Acc
def Model_ACC(model, X_test_first, Y_test_first, BATCH_SIZE):
    y_pred = model.predict(X_test_first, batch_size=BATCH_SIZE)
    y_pred = format_predict(y_pred)
    initial_positive = 0
    predict_right = 0
    predict_wrong = 0

    initial_negtive = 0
    predict_negtive = 0
    predict_right_negtive = 0

    for i in range(len(y_pred)):
        predict_tmp = y_pred[i]
        initial_tmp = Y_test_first[i]

        if int(initial_tmp) == 1:
            initial_positive += 1
        if int(predict_tmp) == 1 and int(initial_tmp) == 1:
            predict_right += 1
        if int(predict_tmp) == 0 and int(initial_tmp) == 1:
            predict_wrong += 1

        if int(initial_tmp) == 0:
            initial_negtive += 1
        if int(predict_tmp) == 1 and int(initial_tmp) == 0:
            predict_negtive += 1
        if int(predict_tmp) == 0 and int(initial_tmp) == 0:
            predict_right_negtive += 1

    if (predict_right + predict_negtive + predict_right_negtive + predict_wrong) != 0:
        res = (predict_right + predict_right_negtive) / (
                predict_right + predict_negtive + predict_right_negtive + predict_wrong)
    else:
        res = 0
    return res


# Sn
def Model_Sn(model, X_test_first, Y_test_first, BATCH_SIZE):
    y_pred = model.predict(X_test_first, batch_size=BATCH_SIZE)
    y_pred = format_predict(y_pred)
    initial_positive = 0
    predict_right = 0
    predict_wrong = 0
    for i in range(len(y_pred)):
        predict_tmp = y_pred[i]
        initial_tmp = Y_test_first[i]

        if int(initial_tmp) == 1:
            initial_positive += 1
        if int(predict_tmp) == 1 and int(initial_tmp) == 1:
            predict_right += 1
        if int(predict_tmp) == 0 and int(initial_tmp) == 1:
            predict_wrong += 1

    if (predict_right + predict_wrong) != 0:
        res = predict_right / (predict_right + predict_wrong)
    else:
        res = 0

    return res


# Sp
def Model_Sp(model, X_test_first, Y_test_first, BATCH_SIZE):
    y_pred = model.predict(X_test_first, batch_size=BATCH_SIZE)
    y_pred = format_predict(y_pred)

    initial_negtive = 0
    predict_negtive = 0
    predict_right_negtive = 0

    for i in range(len(y_pred)):
        predict_tmp = y_pred[i]
        initial_tmp = Y_test_first[i]

        if int(initial_tmp) == 0:
            initial_negtive += 1
        if int(predict_tmp) == 1 and int(initial_tmp) == 0:
            predict_negtive += 1
        if int(predict_tmp) == 0 and int(initial_tmp) == 0:
            predict_right_negtive += 1

    if (predict_right_negtive + predict_negtive) != 0:
        res = predict_right_negtive / (predict_right_negtive + predict_negtive)
    else:
        res = 0

    return res


# MCC
def Model_MCC(model, X_test_first, Y_test_first, BATCH_SIZE):
    y_pred = model.predict(X_test_first, batch_size=BATCH_SIZE)
    y_pred = format_predict(y_pred)
    initial_positive = 0
    predict_right = 0
    predict_wrong = 0

    initial_negtive = 0
    predict_negtive = 0
    predict_right_negtive = 0

    for i in range(len(y_pred)):
        predict_tmp = y_pred[i]
        initial_tmp = Y_test_first[i]

        if int(initial_tmp) == 1:
            initial_positive += 1
        if int(predict_tmp) == 1 and int(initial_tmp) == 1:
            predict_right += 1
        if int(predict_tmp) == 0 and int(initial_tmp) == 1:
            predict_wrong += 1

        if int(initial_tmp) == 0:
            initial_negtive += 1
        if int(predict_tmp) == 1 and int(initial_tmp) == 0:
            predict_negtive += 1
        if int(predict_tmp) == 0 and int(initial_tmp) == 0:
            predict_right_negtive += 1

    fenzi = (predict_right * predict_right_negtive) - (predict_negtive * predict_wrong)
    fenmu = ((predict_right + predict_wrong) * (predict_right + predict_wrong) * (
            predict_right_negtive + predict_negtive) * (predict_right_negtive + predict_wrong)) ** 0.5
    if fenmu != 0:
        res = fenzi / fenmu
    else:
        res = 0

    return res


# AUC
def Model_AUC(model, X_test_first, Y_test_first, BATCH_SIZE):
    Test_y_pred_score = model.predict(X_test_first, batch_size=BATCH_SIZE)
    dataframe = pd.DataFrame({'真实标签': Y_test_first, '预测概率': Test_y_pred_score[:, 0]})
    dataframe.to_csv("/data/lyli/Stomach/stomach-model/Bert_LSTM_DNN_model/34years/stomach_34years_method.csv", index=False, sep=',')
    Test_y_true = np.array(Y_test_first)
    Test_y_scores = np.array(Test_y_pred_score)
    return roc_auc_score(Test_y_true, Test_y_scores)

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

filename1 = "/data/lyli/Stomach/Processed_dataset/enhancer/en_feature_34years_RE_AE_test789.txt"  # 测试集
filename2 = "/data/lyli/Stomach/Processed_dataset/nenhancer/nen_feature_34years_RE_AE_test789.txt"
filename3 = "/data/lyli/Stomach/Processed_dataset/enhancer/window11_test_789.txt"
filename4 = "/data/lyli/Stomach/Processed_dataset/nenhancer/window11_test_789.txt"
X_omics_test, Label12 = load_data(filename1, filename2)
Y_omics_test = []
for i in Label12:
    Y_omics_test.append(i)
Y_omics_test = np.array(Y_omics_test)


sequenceList = read_Fasta(filename3, filename4)
random.seed(75)
random.shuffle(sequenceList)  # shuffle
# print("样本个数",len(sequenceList))

wordVectorSeqList, wordVectorTarList = k_mer_split(sequenceList, k_merList)


# load bert model
config_path = '/data/lyli/Stomach/stomach-model/Bert-embedding/cased_L-12_H-768_A-12/bert_config.json'
checkpoint_path = '/data/lyli/Stomach/stomach-model/Bert-embedding/cased_L-12_H-768_A-12/bert_model.ckpt'
dict_path = '/data/lyli/Stomach/stomach-model/Bert-embedding/cased_L-12_H-768_A-12/vocab.txt'
model_bert = load_trained_model_from_checkpoint(config_path, checkpoint_path)

# Load SeqPoseDict
token_dict = np.load('/data/lyli/Stomach/stomach-model/Bert-embedding/stomach_seq_Dict.npy', allow_pickle=True).item()

tokenizer = OurTokenizer(token_dict)
X_seq_test = [extract_single_features(tokenizer, x['sequence']) for x in sequenceList]
X_seq_test = np.array(X_seq_test)
Y_seq_test = wordVectorTarList




# load model
model_first = Model_TransBERT(X_seq_test.shape[1], 6)
model_first.load_weights(model_weight_filepath)

res = Model_ACC(model_first, [X_seq_test, X_omics_test], Y_seq_test, BATCH_SIZE)
print('The Acc of the first layer model on the test set is:', res)

res = Model_Sn(model_first, [X_seq_test, X_omics_test], Y_seq_test, BATCH_SIZE)
print('The Sn of the first layer model on the test set is:', res)

res = Model_Sp(model_first, [X_seq_test, X_omics_test], Y_seq_test, BATCH_SIZE)
print('The Sp of the first layer model on the test set is:', res)

res = Model_MCC(model_first, [X_seq_test, X_omics_test], Y_seq_test, BATCH_SIZE)
print('The MCC of the first layer model on the test set is:', res)

res = Model_AUC(model_first, [X_seq_test, X_omics_test], Y_seq_test, BATCH_SIZE)
print('The AUC of the first layer model on the test set is:', res)






