
import numpy as np
from keras.models import *
from keras.layers import Dense, Dropout, LSTM, GRU,Embedding, Bidirectional, Concatenate
from keras.layers import Input, Dense, Dropout, BatchNormalization, LSTM, Bidirectional
import os
from keras import backend as K
import random
from keras.callbacks import ModelCheckpoint
from keras.optimizers import Adam

# Hyperparameter settings
from sklearn.model_selection import StratifiedKFold

BATCH_SIZE = 64    # best batch-size
EPCOHS = 10        # best epoch
RNN_HIDDEN_DIM = 64  # Number of neural units in the lstm layer
DROPOUT_RATIO_1 = 0.5  # proportion of neurones not used for training
DROPOUT_RATIO_2 = 0.5  # proportion of neurones not used for training
DROPOUT_RATIO_3 = 0.1


checkpoint_dir = '../Liver/liver-model/model/25years/final_model_25years'  # Store model files
if not os.path.exists(checkpoint_dir):  # if the path does not exist, generate it
    os.mkdir(checkpoint_dir)

model_weight_filepath = checkpoint_dir + "/liver_25years_model.h5"  # First layer model
model_structure_filepath = checkpoint_dir + "/liver_25years_model.json"  # First layer model


# load multi_omics data features
def load_omics_feature(filename1, filename2):
    filelist = [filename1, filename2]
    dataMat = []
    for file in filelist:
        # file0 = open(file)
        file0 = np.loadtxt(file, encoding='utf_8_sig', delimiter=' ')
        # for line in file0.readlines():
        for line in file0:
            # curLine = line.strip().split(" ")
            curLine = list(line)
            # print(curLine)
            floatLine = map(float, curLine)  # 这里使用的是map函数直接把数据转化成为float类型
            floatLine = list(floatLine)
            dataMat.append(floatLine[0:7])

    # print(np.shape(dataMat))

    random.seed(75)
    random.shuffle(dataMat)  # shuffle

    labelMat = []
    for i in dataMat:
        labelMat.append(int(i[-1]))
        del (i[-1])

    # scaler = MinMaxScaler()  # 线性归一化转换
    # dataMat = scaler.fit_transform(dataMat)  # 转换数据集

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
        # attention_mul_bilstm = attention_3d_block(bn, input_length)
        drop_1 = Dropout(DROPOUT_RATIO_1)(bn)
        bi_lstm_2 = Bidirectional(LSTM(rnn_hidden_dim))(drop_1)
        drop_2 = Dropout(DROPOUT_RATIO_2)(bi_lstm_2)
        seq_dense = Dense(1, activation='sigmoid')(drop_2)
        print("drop_2格式")
        print(np.shape(drop_2))
        # (?, 128)
        # (?, 1)

        omics_dense1 = Dense(50, name="chip_first_layer")(inputs_omics)
        omics_dense2 = Dropout(DROPOUT_RATIO_3)(omics_dense1)
        omics_dense3 = Dense(4, name="chip_second_layer", activation="relu")(omics_dense2)
        omics_dense4 = BatchNormalization()(omics_dense3)
        # omics_dense1 = Dense(128, activation='relu')(inputs_omics)
        # omics_dense2 = Dense(64, activation='relu')(omics_dense1)
        # omics_dense3 = Dense(64, activation='relu')(omics_dense2)
        # omics_dense4 = Dense(128, activation='relu')(omics_dense3)

        print(np.shape(omics_dense4))
        print("omics_dense4格式")
        # x = merge([drop_2, omics_dense4], mode='concat', concat_axis=2)
        x = Concatenate(axis=1)([seq_dense, omics_dense4])
        # x = Flatten()(x)
        # dense1_ = Dense(512, activation='relu')
        # dense1 = dense1_(x)
        # dense2 = Dense(256, activation='relu')(dense1)
        # x = Dropout(0.5)(dense2)
        # dense3 = Dense(128, activation='relu')(x)
        # pred_output = Dense(1, activation='sigmoid')(dense3)

        pred_output = Dense(1, activation='sigmoid')(x)
        model = Model(input=[inputs_seq1, inputs_omics], output=[pred_output])

        adam = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001 / EPCOHS)
        model.compile(loss='binary_crossentropy', optimizer=adam, metrics=["accuracy"])

        # dense = Dense(1, activation='sigmoid')(drop_2)
        # model = Model(inputs_seq1, inputs_omics, dense)
        # model.compile('adam', 'binary_crossentropy', metrics=['accuracy'])

        return model

    return create_lstm(input_seq_length, input_omics_length)


# load sequence's word vector
X_seq_train = np.loadtxt('../Liver/liver-model/model/Bert-embedding/liver_traindataset_seq_wordVector.txt', delimiter=',').reshape((5074, 512, 768))
X_seq_train = np.array(X_seq_train)
# Y_train = wordVectorTarList

# load omics feature
# training dataset
filename1 = "../Liver/Processed_dataset/enhancer/en_feature_25years_RE_AE_train2537.txt"
filename2 = "../Liver/Processed_dataset/nenhancer/nen_feature_25years_train2537.txt"
X_omics_train, Y_omics_train = load_omics_feature(filename1, filename2)
Y_omics_train = np.array(Y_omics_train)



# training Model
if __name__ == '__main__':
    model = Model_TransBERT(X_seq_train.shape[1], 6)  # Get the first layer model
    # saveBestModel = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, mode='max')
    # callbacks_list = [saveBestModel]
    print('Fitting first model...')
    # history = model.fit([X_seq_train, X_omics_train], Y_train, batch_size=BATCH_SIZE, epochs=EPCOHS, callbacks=callbacks_list, vvalidation_data=[[X_seq_valid, X_omics_valid], Y_valid], verbose=1)  # Start training
    model.fit([X_seq_train, X_omics_train], Y_omics_train, batch_size=BATCH_SIZE, epochs=EPCOHS, verbose=1)  # Start training

    # serialize model to JSON
    model_json = model.to_json()
    with open(model_structure_filepath, "w") as json_file:
        json_file.write(model_json)
        # serialize weights to HDF5
    model.save_weights(model_weight_filepath)
    print("Saved model to disk")


