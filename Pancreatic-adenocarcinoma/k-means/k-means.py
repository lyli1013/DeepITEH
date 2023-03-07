



import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt

def cluster_RE_AE(filename1, filename2, filename3):
    enhancer = np.loadtxt(filename1, encoding='utf_8_sig', delimiter=' ')
    nenhancer = np.loadtxt(filename2, encoding='utf_8_sig', delimiter=' ')
    X = np.concatenate((enhancer, nenhancer), axis=0)

    min_max_scale = MinMaxScaler()
    X = min_max_scale.fit_transform(X)
    kmeans = KMeans(n_clusters = 2, init = 'k-means++', random_state = 5).fit(X)
    colors = ['red', 'green', 'blue', 'black', 'yellow', 'purple', 'grey', 'pink', 'orange', 'brown', 'navy', 'cyan']

    class1 = []
    class2 = []
    for i, cluster in enumerate(kmeans.labels_):
        plt.scatter(X[i][0], X[i][1], color=colors[cluster])
        # print(X[i], colors[cluster], cluster, kmeans.labels_)
        if colors[cluster] == 'red':
            class1.append(i)
            # print(X[i])
        if colors[cluster] == 'green':
            class2.append(i)
            # print(X[i])
    # 哪个类别中数量少哪个是RE
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

    return 0


lung_enhancer = r"..\Pancreatic-adenocarcinoma\PAAD-Multi-omics_data\Panc1\Enhancer_feature\en_feature_835_Panc1.txt"
lung_nenhancer = r"..\Pancreatic-adenocarcinoma\NoEnhancer\nen_feature_835_Panc1.txt"
Liver_16years_RE_index_file = r"..\Pancreatic-adenocarcinoma\Processed_dataset-new\Panc1_80_index.txt"
cluster_RE_AE(lung_enhancer, lung_nenhancer, Liver_16years_RE_index_file)



# 在每个enhancer/nenhancer后添加RE(1)或RE/NE(0)标签
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
        if n < 1178:
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

# 按照投票后的RE，给从LIHC的组学数据中提取的增强子和非增强子组学特征添加RE_AE_lable
Liver_31years_enhancer_omics_feature_lable_file = r"..\Pancreatic-adenocarcinoma\Processed_dataset-new\enhancer\en_feature_835_Panc1_RE_AE.txt"
Liver_31years_nenhancer_omics_feature_lable_file = r"..\Pancreatic-adenocarcinoma\Processed_dataset-new\nenhancer\nen_feature_835_Panc1_RE_AE.txt"
RE_AE_lable(lung_enhancer, lung_nenhancer, Liver_16years_RE_index_file, Liver_31years_enhancer_omics_feature_lable_file, Liver_31years_nenhancer_omics_feature_lable_file)
