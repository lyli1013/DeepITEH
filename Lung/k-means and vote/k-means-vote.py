'''
# 在从每套组学数据中提取增强子和非增强子的组学信号后，将增强子和非增强子放在一起按照每套组学特征（）进行聚类，聚成两类，即RE、AE/NE
# 有几套组学数据就是聚出几套RE
# 然后将几套的RE进行投票，得到投票后的RE
# 用投票后得到的RE和AE/NE，对每套增强子和非增强子打RE/AE的标签，即在每套增强子和非增强子组学特征矩阵中添加1/0
# 再在每套增强子和非增强子的组学特征中打增强子和非增强子的标签
'''


import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt

# 对每套增强子和非增强子的组学特征聚类，得到每套数据的RE下标
def cluster_RE_AE(filename1, filename2, filename3):
    enhancer = np.loadtxt(filename1, encoding='utf_8_sig', delimiter=' ')
    nenhancer = np.loadtxt(filename2, encoding='utf_8_sig', delimiter=' ')
    X = np.concatenate((enhancer, nenhancer), axis=0)
    min_max_scale = MinMaxScaler()
    X = min_max_scale.fit_transform(X)
    kmeans = KMeans(n_clusters=2, init='k-means++', random_state=5).fit(X)
    colors = ['red', 'green', 'blue', 'black', 'yellow', 'purple', 'grey', 'pink', 'orange', 'brown', 'navy', 'cyan']
    class1 = []
    class2 = []
    # class3 = []
    # class4 = []
    # class5 = []
    # class6 = []
    # class7 = []
    # class8 = []
    # class9 = []
    # class10 = []
    # class11 = []
    # class12 = []
    for i, cluster in enumerate(kmeans.labels_):
        plt.scatter(X[i][0], X[i][1], color=colors[cluster])
        # print(X[i], colors[cluster], cluster, kmeans.labels_)
        if colors[cluster] == 'red':
            class1.append(i)
            # print(X[i])
        if colors[cluster] == 'green':
            class2.append(i)
            # print(X[i])
        # if colors[cluster] == 'blue':
        #     class3.append(i)
        #     # print(X[i])
        # if colors[cluster] == 'black':
        #     class4.append(i)
        #     # print(X[i])
        # if colors[cluster] == 'yellow':
        #     class5.append(i)
        #     # print(X[i])
        # if colors[cluster] == 'purple':
        #     class6.append(i)
        #     # print(X[i])
        # if colors[cluster] == 'grey':
        #     class7.append(i)
        #     # print(X[i])
        # if colors[cluster] == 'pink':
        #     class8.append(i)
        #     # print(X[i])
        # if colors[cluster] == 'orange':
        #     class9.append(i)
        #     # print(X[i])
        # if colors[cluster] == 'brown':
        #     class1.append(i)
        #     # print(X[i])
        # if colors[cluster] == 'navy':
        #     class1.append(i)
        #     # print(X[i])
        # if colors[cluster] == 'cyan':
        #     class1.append(i)
        #     # print(X[i])
    print(len(class1))
    print(len(class2))
    # print(len(class3))
    # print(len(class4))
    # print(len(class5))
    # print(len(class6))
    # print(len(class7))
    # print(len(class8))
    # print(len(class9))
    # print(len(class10))
    # print(len(class11))
    # print(len(class12))
    # # print(class2)

    file = open(filename3, 'w')
    if len(class1) < len(class2):
        RE_nums = len(class1)
        AE_nums = len(class2)
        print("RE:", RE_nums)
        print("AE:", AE_nums)
        for j in class1:
            file.write(str(j) + '\n')
        file.close()
    else:
        RE_nums = len(class2)
        AE_nums = len(class1)
        print("RE:", RE_nums)
        print("AE:", AE_nums)
        for j in class2:
            file.write(str(j) + '\n')
        file.close()
    # plt.show()

    return 0

# 五个RE文件投票，得到投票后的RE下标文件
def RE_vote(filename1, filename2, filename3, filename4, filename5, filename6):
    years_3 = []
    with open(filename1, 'r')as f:
        for line in f:
            years_3.append(int(line.strip()))    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开
    # print(years_16)
    years_37 = []
    with open(filename2, 'r')as f:
        for line in f:
            years_37.append(int(line.strip()))    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

    years_51 = []
    with open(filename3, 'r')as f:
        for line in f:
            years_51.append(int(line.strip()))    # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

    years_54 = []
    with open(filename4, 'r') as f:
        for line in f:
            years_54.append(int(line.strip()))  # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

    years_59 = []
    with open(filename5, 'r') as f:
        for line in f:
            years_59.append(int(line.strip()))  # strip()表示删除掉数据中的换行符，split（‘ ’）则是数据中遇到空格就隔开

    res3_37_51 = list(set(years_3) & set(years_37) & set(years_51))
    res3_37_54 = list(set(years_3) & set(years_37) & set(years_54))
    res3_37_59 = list(set(years_3) & set(years_37) & set(years_59))
    res3_51_54 = list(set(years_3) & set(years_51) & set(years_54))
    res3_51_59 = list(set(years_3) & set(years_51) & set(years_59))
    res3_54_59 = list(set(years_3) & set(years_54) & set(years_59))
    res37_51_54 = list(set(years_37) & set(years_51) & set(years_54))
    res37_54_59 = list(set(years_37) & set(years_54) & set(years_59))
    res37_51_59 = list(set(years_37) & set(years_51) & set(years_59))
    res51_54_59 = list(set(years_51) & set(years_54) & set(years_59))
    print(len(res3_37_51))
    print(len(res3_37_54))
    print(len(res3_37_59))
    print(len(res3_51_54))
    print(len(res3_51_59))
    print(len(res3_54_59))
    print(len(res37_51_54))
    print(len(res37_54_59))
    print(len(res37_51_59))
    print(len(res51_54_59))

    res3_37_51_54 = list(set(years_3) & set(years_37) & set(years_51) & set(years_54))
    res3_37_51_59 = list(set(years_3) & set(years_37) & set(years_51) & set(years_59))
    res3_37_54_59 = list(set(years_3) & set(years_37) & set(years_54) & set(years_59))
    res3_51_54_59 = list(set(years_3) & set(years_51) & set(years_54) & set(years_59))
    res37_51_54_59 = list(set(years_37) & set(years_51) & set(years_54) & set(years_59))
    res3_37_51_54_59 = list(set(years_3) & set(years_37) & set(years_51) & set(years_54) & set(years_59))
    print(len(res3_37_51_54))
    print(len(res3_37_51_59))
    print(len(res3_37_54_59))
    print(len(res3_51_54_59))
    print(len(res37_51_54_59))
    print(len(res3_37_51_54_59))



    ret1=list(set(res3_37_51).union(set(res3_37_54)).union(set(res3_37_59)).union(set(res3_51_54)).union(set(res3_51_59))
              .union(set(res3_54_59)).union(set(res37_51_54)).union(set(res37_54_59)).union(set(res37_51_59)).union(set(res51_54_59))
              .union(set(res3_37_51_54)).union(set(res3_37_51_59)).union(set(res3_37_54_59)).union(set(res3_51_54_59)).union(set(res37_51_54_59)).union(set(res3_37_51_54_59)))
    ret1.sort()
    print(ret1)
    print(len(ret1))

    # res16_25_31 = list(set(years_3) & set(years_37) & set(years_51))
    # res16_25 = list(set(years_3) & set(years_37))
    # res16_31 = list(set(years_3) & set(years_51))
    # res25_31 = list(set(years_37) & set(years_51))
    # print(len(res16_25))
    # print(len(res16_31))
    # print(len(res25_31))
    # print(len(res16_25_31))
    # ret1 = list(set(res16_25_31).union(set(res16_25)).union(set(res16_31)).union(set(res25_31)))
    # ret1.sort()
    # print(ret1)
    # print(len(ret1))
    file = open(filename6, 'w')
    for j in ret1:
        file.write(str(j) + '\n')
    file.close()

    return 0

# 用投票后的RE下标文件分别给三套/几套的增强子和非增强子组学特征进行标注1/0，并添加增强子和非增强子的本身标签
def RE_AE_lable(filename1, filename2, filename3, filename4, filename5):
    enhancer = np.loadtxt(filename1, encoding='utf_8_sig', delimiter=' ')
    nenhancer = np.loadtxt(filename2, encoding='utf_8_sig', delimiter=' ')
    enhancer_nums = len(enhancer)
    nenhancer_nums = len(nenhancer)
    X = np.concatenate((enhancer, nenhancer), axis=0)

    X_new = []
    index = []
    file_index = open(filename3, 'r')
    for i in file_index:
        la = int(i)
        index.append(la)
    RE_nums = len(index)
    n = 0
    m = 0
    for enhancer in X:
        if m < RE_nums and n == int(index[m]):
            enhancer_list = list(enhancer)
            enhancer_list.append(1)
            m = m + 1
        else:
            enhancer_list = list(enhancer)
            enhancer_list.append(0)
        if n < 6481:
            enhancer_list.append(1)
        else:
            enhancer_list.append(0)
        n = n + 1
        X_new.append(enhancer_list)

    X_new = np.array(X_new)
    enhancer_RE_AE = []
    nenhancer_RE_AE = []
    for sample in range(enhancer_nums):
        enhancer_RE_AE.append(X_new[sample])
    enhancer_RE_AE = np.array(enhancer_RE_AE)

    for sample in range(enhancer_nums, (enhancer_nums + nenhancer_nums)):
        nenhancer_RE_AE.append(X_new[sample])
    nenhancer_RE_AE = np.array(nenhancer_RE_AE)

    np.savetxt(filename4, enhancer_RE_AE, encoding='utf_8_sig', fmt='%.3f', delimiter=' ')
    np.savetxt(filename5, nenhancer_RE_AE, encoding='utf_8_sig', fmt='%.3f', delimiter=' ')

    return 0

# 1、按照从Liver_16years的组学数据中提取的enhancer和nenhancer生物信号聚类，获得RE
Lung_3years_enhancer_omics_feature_file = r"..\Lung\Lung-Multi-omics_data\lung_3years\Enhancer_feature\en_feature_6481_3years.txt"
Lung_3years_nenhancer_omics_feature_file = r"..\NoEnhancer\NoEnhancer_3years_feature_6481_6001.txt"
Lung_3years_RE_index_file = r"..\Lung\Processed_dataset\feature_3years_index.txt"
cluster_RE_AE(Lung_3years_enhancer_omics_feature_file, Lung_3years_nenhancer_omics_feature_file, Lung_3years_RE_index_file)

# 按照从Liver_25years的组学数据中提取的enhancer和nenhancer生物信号聚类，获得RE
Lung_37years_enhancer_omics_feature_file = r"..\Lung\Lung-Multi-omics_data\upper_lobe_of_left_lung_37years\Enhancer_feature\en_feature_6481_37years.txt"
Lung_37years_nenhancer_omics_feature_file = r"..\NoEnhancer\NoEnhancer_37years_feature_6481_6001.txt"
Lung_37years_RE_index_file = r"..\Lung\Processed_dataset\feature_37years_index.txt"
cluster_RE_AE(Lung_37years_enhancer_omics_feature_file, Lung_37years_nenhancer_omics_feature_file, Lung_37years_RE_index_file)

# 按照从Liver_31years的组学数据中提取的enhancer和nenhancer生物信号聚类，获得RE
Lung_51years_enhancer_omics_feature_file = r"..\Lung\Lung-Multi-omics_data\upper_lobe_of_left_lung_51years\Enhancer_feature\en_feature_6481_51years.txt"
Lung_51years_nenhancer_omics_feature_file = r"..\NoEnhancer\NoEnhancer_51years_feature_6481_6001.txt"
Lung_51years_RE_index_file = r"..\Lung\Processed_dataset\feature_51years_index.txt"
cluster_RE_AE(Lung_51years_enhancer_omics_feature_file, Lung_51years_nenhancer_omics_feature_file, Lung_51years_RE_index_file)

# 按照从Liver_31years的组学数据中提取的enhancer和nenhancer生物信号聚类，获得RE
Lung_54years_enhancer_omics_feature_file = r"..\Lung\Lung-Multi-omics_data\upper_lobe_of_left_lung_54years\Enhancer_feature\en_feature_6481_54years.txt"
Lung_54years_nenhancer_omics_feature_file = r"..\NoEnhancer\NoEnhancer_54years_feature_6481_6001.txt"
Lung_54years_RE_index_file = r"..\Lung\Processed_dataset\feature_54years_index.txt"
cluster_RE_AE(Lung_54years_enhancer_omics_feature_file, Lung_54years_nenhancer_omics_feature_file, Lung_54years_RE_index_file)

# 按照从Liver_31years的组学数据中提取的enhancer和nenhancer生物信号聚类，获得RE
Lung_59years_enhancer_omics_feature_file = r"..\Lung\Lung-Multi-omics_data\lower_lobe_of_left_lung_59years_hg38\hg19\Enhancer_feature\en_feature_6481_59years.txt"
Lung_59years_nenhancer_omics_feature_file = r"..\NoEnhancer\NoEnhancer_59years_feature_6481_6001.txt"
Lung_59years_RE_index_file = r"..\Lung\Processed_dataset\feature_59years_index.txt"
cluster_RE_AE(Lung_59years_enhancer_omics_feature_file, Lung_59years_nenhancer_omics_feature_file, Lung_59years_RE_index_file)

# 2、将以上三种RE进行vote，获得真实性更高的RE
Lung_RE_vote_index_file = r"..\Lung\Processed_dataset\lung_feature_RE_index.txt"
RE_vote(Lung_3years_RE_index_file, Lung_37years_RE_index_file, Lung_51years_RE_index_file, Lung_54years_RE_index_file, Lung_59years_RE_index_file, Lung_RE_vote_index_file)

# 3、按照投票后的RE，给从16years—liver的组学数据中提取的增强子和非增强子组学特征添加RE_AE_lable
Lung_3years_enhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\enhancer\en_feature_6481_3years_RE_AE.txt"
Lung_3years_nenhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\nenhancer\NoEnhancer_3years_feature_6481_6001_RE_AE.txt"
RE_AE_lable(Lung_3years_enhancer_omics_feature_file, Lung_3years_nenhancer_omics_feature_file, Lung_RE_vote_index_file, Lung_3years_enhancer_omics_feature_lable_file, Lung_3years_nenhancer_omics_feature_lable_file)


Lung_37years_enhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\enhancer\en_feature_6481_37years_RE_AE.txt"
Lung_37years_nenhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\nenhancer\NoEnhancer_37years_feature_6481_6001_RE_AE.txt"
RE_AE_lable(Lung_37years_enhancer_omics_feature_file, Lung_37years_nenhancer_omics_feature_file, Lung_RE_vote_index_file, Lung_37years_enhancer_omics_feature_lable_file, Lung_37years_nenhancer_omics_feature_lable_file)

Lung_51years_enhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\enhancer\en_feature_6481_51years_RE_AE.txt"
Lung_51years_nenhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\nenhancer\NoEnhancer_51years_feature_6481_6001_RE_AE.txt"
RE_AE_lable(Lung_51years_enhancer_omics_feature_file, Lung_51years_nenhancer_omics_feature_file, Lung_RE_vote_index_file, Lung_51years_enhancer_omics_feature_lable_file, Lung_51years_nenhancer_omics_feature_lable_file)

Lung_54years_enhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\enhancer\en_feature_6481_54years_RE_AE.txt"
Lung_54years_nenhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\nenhancer\NoEnhancer_54years_feature_6481_6001_RE_AE.txt"
RE_AE_lable(Lung_54years_enhancer_omics_feature_file, Lung_54years_nenhancer_omics_feature_file, Lung_RE_vote_index_file, Lung_54years_enhancer_omics_feature_lable_file, Lung_54years_nenhancer_omics_feature_lable_file)

Lung_59years_enhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\enhancer\en_feature_6481_59years_RE_AE.txt"
Lung_59years_nenhancer_omics_feature_lable_file = r"..\Lung\Processed_dataset\nenhancer\NoEnhancer_59years_feature_6481_6001_RE_AE.txt"
RE_AE_lable(Lung_59years_enhancer_omics_feature_file, Lung_59years_nenhancer_omics_feature_file, Lung_RE_vote_index_file, Lung_59years_enhancer_omics_feature_lable_file, Lung_59years_nenhancer_omics_feature_lable_file)
