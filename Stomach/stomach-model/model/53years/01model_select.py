import numpy as np
from keras.models import *
from keras.layers import Input, Dense, Dropout, BatchNormalization, LSTM, Bidirectional, Concatenate
from keras import backend as K
import random
from keras.optimizers import Adam
from sklearn.model_selection import StratifiedKFold

# Hyperparameter settings
RNN_HIDDEN_DIM = 64  # Number of neural units in the lstm layer
DROPOUT_RATIO_1 = 0.5  # proportion of neurones not used for training
DROPOUT_RATIO_2 = 0.5  # proportion of neurones not used for training
DROPOUT_RATIO_3 = 0.1


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


# load train_val sequence's word vector
X_seq_train_val = np.loadtxt('/data/lyli/Stomach/stomach-model/Bert-embedding/stomach_seq_train_6308_wordVector.txt', delimiter=',').reshape((6308, 512, 768))
X_seq_train_val = np.array(X_seq_train_val)

# load omics feature
filename1 = "/data/lyli/Stomach/Processed_dataset/enhancer/en_feature_53years_RE_AE_train3154.txt"
filename2 = "/data/lyli/Stomach/Processed_dataset/nenhancer/nen_feature_53years_RE_AE_train3154.txt"
X_omics_train_val, Label = load_omics_feature(filename1, filename2)
Label = np.array(Label)


# Select model/hyperparameter
best_ACC = 0.0
folds = StratifiedKFold(n_splits=5, shuffle=False)
count_iter = 1  # Record the number of cycles
for BATCH_SIZE in [16,32,64]:
# for BATCH_SIZE in [16]:
# for BATCH_SIZE in [32]:
# for BATCH_SIZE in [64]:
    for EPCOHS in [3, 4, 5, 6, 7, 8, 9, 10]:
        fold_nums = 1
        ACC = []
        for train_index, valid_index in folds.split(X_omics_train_val, Label):
            # train_index,test_index是数组索引，如[1 2 3]
            print("This is Epoch="+str(EPCOHS)+";"+"Batch_size="+str(BATCH_SIZE))
            print("This fold_nums is:", fold_nums)
            K.clear_session()  # Release graph
            X_seq_train = X_seq_train_val[train_index]  # X_train_first是数组
            X_omics_train = X_omics_train_val[train_index]
            Y_train = Label[train_index]

            X_seq_valid = X_seq_train_val[valid_index]  # X_train_first是数组
            X_omics_valid = X_omics_train_val[valid_index]
            Y_valid = Label[valid_index]

            X_seq_train = np.array(X_seq_train)
            X_omics_train = np.array(X_omics_train)
            X_seq_valid = np.array(X_seq_valid)
            X_omics_valid = np.array(X_omics_valid)

            # select model
            model = Model_TransBERT(X_seq_train.shape[1], 6)  # Get the first layer model
            # saveBestModel = ModelCheckpoint(filepath, monitor='val_acc', verbose=1, save_best_only=True, mode='max')
            # callbacks_list = [saveBestModel]
            print('Fitting first model...')
            # history = model.fit([X_seq_train, X_omics_train], Y_train, batch_size=BATCH_SIZE, epochs=EPCOHS, callbacks=callbacks_list, vvalidation_data=[[X_seq_valid, X_omics_valid], Y_valid], verbose=1)  # Start training
            model.fit([X_seq_train, X_omics_train], Y_train, batch_size=BATCH_SIZE, epochs=EPCOHS, verbose=1)  # Start training
            loss, accuracy = model.evaluate([X_seq_valid, X_omics_valid], Y_valid, batch_size=BATCH_SIZE)
            print('\ntest loss', loss)
            print('accuracy', accuracy)
            ACC.append(accuracy)

            K.clear_session()
            del X_seq_train
            del X_omics_train
            del Y_train
            del X_seq_valid
            del X_omics_valid
            del Y_valid

            count_iter += 1
            fold_nums += 1

        ACC1 = np.array(ACC)
        ACC2 = ACC1.mean()
        print("acc", ACC2)
        if ACC2 > best_ACC:
            best_ACC = ACC2
            best_parameter = {'batch_size': BATCH_SIZE, 'epochs': EPCOHS}
print("Best score:", best_ACC)
print("best_parameter", best_parameter)
