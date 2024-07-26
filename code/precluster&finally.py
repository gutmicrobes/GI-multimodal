root = r'E:\\BaiduNetdiskDownload\\王瑜多_毕业论文代码数据\\dyang\\pre-data'

import os
os.chdir(root)
import pandas as pd
mirna = pd.read_csv('dis_mirna0627.csv')
mrna = pd.read_csv('dis_mrna0627.csv')
image = pd.read_csv('dis_image0627.csv')
immune=pd.read_csv('dis_immune0627.csv')

mrna = mrna.iloc[:,1:]
mirna = mirna.iloc[:,1:]
image = image.iloc[:,1:]
immune=immune.iloc[:,1:]

from sklearn_extra.cluster import KMedoids
import numpy as np
from sklearn.metrics import silhouette_score
X=mirna.values
X=mrna.values
X=image.values
X=immune.values

# 计算不同K值的轮廓系数
silhouette_scores = []
K_max = 10  # 最大聚类数
for k in range(2, K_max+1):
    kmedoids = KMedoids(n_clusters=k, random_state=0).fit(X)
    labels = kmedoids.labels_
    silhouette_scores.append(silhouette_score(X, labels))

# 找到轮廓系数最高的K值
best_k = range(2, K_max+1)[silhouette_scores.index(max(silhouette_scores))]
print("Best K by silhouette score:", best_k)
#mirna 2  mrna 2

import matplotlib.pyplot as plt

# 计算总内部距离
total_intra_distances = []
for k in range(2, K_max+1):
    kmedoids = KMedoids(n_clusters=k, random_state=0).fit(X)
    total_intra_distances.append(kmedoids.inertia_)

# 绘制肘部图
plt.figure(figsize=(8, 4))
plt.plot(range(2, K_max+1), total_intra_distances, 'bo-')
plt.xlabel('Number of clusters K')
plt.ylabel('Total intra-cluster distance')
plt.title('Elbow Method For Optimal K')
plt.show()

#mirna 5  mrna 5 image 4 immune 3
kmedoids = KMedoids(n_clusters=5, random_state=0)
mirna_labels = kmedoids.fit_predict(mirna.values)

mrna_labels = kmedoids.fit_predict(mrna.values)

kmedoids = KMedoids(n_clusters=4, random_state=0)
image_labels = kmedoids.fit_predict(image.values)

kmedoids = KMedoids(n_clusters=3, random_state=0)
immune_labels = kmedoids.fit_predict(immune.values)

# 保存结果
np.savetxt('mirna_labels.txt', mirna_labels, fmt='%d')
np.savetxt('mrna_labels.txt', mrna_labels, fmt='%d')
np.savetxt('image_labels.txt', image_labels, fmt='%d')
np.savetxt('immune_labels.txt', immune_labels, fmt='%d')

import math
def power_mat(data):
    mean_dis = []  #定义样本与其他样本的平均距离的列表
    n = data.shape[0]
    for i in range(n):
        sum1 = 0
        for j in range(n):
            if i!=j:
                sum1+=data.iloc[i,j]
        mean_dis.append(sum1/(n-1))
    p = pd.DataFrame(np.zeros([n,n]))
    u = 0.5
    for i in range(n):
        for j in range(n):
            v = (mean_dis[i]+mean_dis[j]+data.iloc[i,j])/3
            p.iloc[i,j] = math.exp(-pow(data.iloc[i,j],2)/(u*v))
    p0 = pd.DataFrame(np.zeros([n,n]))
    for i in range(n):
        sum2 = sum(p.iloc[i,:])
        for j in range(n):
            if i == j:
                p0.iloc[i,j] = 0.5
            else:
                p0.iloc[i,j] = p.iloc[i,j]/(2*(sum2-p.iloc[i,i]))
    return p,p0
#p0是权重矩阵，p是高斯核矩阵


image_center,image_power = power_mat(image)
mrna_center,mrna_power = power_mat(mrna)
mirna_center,mirna_power = power_mat(mirna)
immune_center,immune_power = power_mat(immune)

#写入文件
image_center.to_csv('image_center.csv',index=False)
image_power.to_csv('image_power.csv',index=False)
mrna_center.to_csv('mrna_center.csv',index=False)
mrna_power.to_csv('mrna_power.csv',index=False)
mirna_center.to_csv('mirna_center.csv',index=False)
mirna_power.to_csv('mirna_power.csv',index=False)
immune_center.to_csv('immune_center.csv',index=False)
immune_power.to_csv('immune_power.csv',index=False)



def core_mat(data, data_list):
    n = data.shape[0]
    S = pd.DataFrame(np.zeros([n, n]))
    for i in range(n):
        sum2 = 0
        for j in range(n):
            if data_list[i] == data_list[j]:
                sum2 += data.iloc[i, j]
        for j in range(n):
            if data_list[i] == data_list[j]:
                S.iloc[i, j] = data.iloc[i, j] / sum2

    return S

mrna_core = core_mat(mrna_center,mrna_labels)
mirna_core = core_mat(mirna_center,mirna_labels)
image_core = core_mat(image_center,image_labels)
immune_core = core_mat(immune_center,immune_labels)

mrna_core.to_csv('mrna_core.csv',index=False)
mirna_core.to_csv('mirna_core.csv',index=False)
image_core.to_csv('image_core.csv',index=False)
immune_core.to_csv('immune_core.csv',index=False)


#标准化
def custom_normalize(W):
    n = W.shape[0]  # 矩阵的大小
    P = np.zeros_like(W)  # 初始化P矩阵

    for i in range(n):
        row_sum = np.sum(W[i, :]) - W[i, i]  # 计算除对角线外的行和
        P[i, i] = 0.5  # 设置对角线为0.5
        for j in range(n):
            if i != j:
                P[i, j] = W[i, j] / (2 * row_sum)  # 设置非对角线元素

    return P

def snf_update(P_list, S_list):
    """ Update each P matrix using all other P matrices, weighted by corresponding S matrices. """
    new_P_list = []
    m = len(P_list)  # Number of data types
    for i in range(m):
        sum_P = np.zeros_like(P_list[i])
        for j in range(m):
            if i != j:
                sum_P += S_list[j] @ P_list[j] @ S_list[j].T
        new_P_list.append(custom_normalize(sum_P / (m - 1)))  # Normalize and average the sum
    return new_P_list

def snf(P_list, S_list, iterations=20):
    """ Perform SNF for a given number of iterations and fuse the matrices. """
    for _ in range(iterations):
        P_list = snf_update(P_list, S_list)
    # Fusion of all P matrices
    fused_P = np.sum(P_list, axis=0) / len(P_list)
    return custom_normalize(fused_P)

def snf(P_list, S_list, iterations=20):
    """ Perform SNF, ensure inputs are NumPy arrays """
    convergence_diffs = []
    for iteration in range(iterations):
        new_P_list = snf_update(P_list, S_list)
        diff = max(np.linalg.norm(P_list[i] - new_P_list[i], 'fro') for i in range(len(P_list)))
        convergence_diffs.append(diff)
        P_list = new_P_list
    fused_P = np.sum(P_list, axis=0) / len(P_list)
    return custom_normalize(fused_P), convergence_diffs

#查看imae_core的是不是dataframe
print(type(image_core))

def snf_update(P_list, S_list):
    """Update each P matrix using all other P matrices, weighted by corresponding S matrices, and return differences."""
    new_P_list = []
    diffs = []
    m = len(P_list)  # Number of data types
    for i in range(m):
        sum_P = np.zeros_like(P_list[i])
        for j in range(m):
            if i != j:
                sum_P += S_list[j] @ P_list[j] @ S_list[j].T

        normalized_P = custom_normalize(sum_P / (m - 1))
        diff = np.linalg.norm(P_list[i] - normalized_P, 'fro')
        diffs.append(diff)
        new_P_list.append(normalized_P)
    return new_P_list, diffs


def snf(P_list, S_list, iterations=4):
    """Perform SNF for a given number of iterations and fuse the matrices, tracking individual diffs."""
    convergence_diffs = [[] for _ in range(len(P_list))]  # List of lists to track diffs for each P matrix

    for _ in range(iterations):
        new_P_list, diffs = snf_update(P_list, S_list)
        for idx, diff in enumerate(diffs):
            convergence_diffs[idx].append(diff)
        P_list = new_P_list

    fused_P = np.sum(P_list, axis=0) / len(P_list)
    fused_P = custom_normalize(fused_P)

    return fused_P, convergence_diffs


def snf(P_list, S_list, tolerance=1e-4, max_iterations=50):
    """Perform SNF until all matrices converge below a tolerance or reach max iterations."""
    convergence_diffs = []  # Initialize list to store differences for convergence check
    iteration = 0
    converged = False

    while not converged and iteration < max_iterations:
        new_P_list, diffs = snf_update(P_list, S_list)
        P_list = new_P_list
        iteration += 1
        convergence_diffs.append(diffs)  # Store differences of each iteration

        # Check if all diffs are below the tolerance
        if all(diff < tolerance for diff in diffs):
            converged = True

    fused_P = np.sum(P_list, axis=0) / len(P_list)
    return fused_P, convergence_diffs  # Make sure to return both values


# Assuming P_list and S_list are initialized and appropriate for the SNF function

# 把image_core,mrna_core,mirna_core,immune_core合并成一个list
p_list = [image_power.values,mrna_power.values,mirna_power.values,immune_power.values]
S_list = [image_core.values,mrna_core.values,mirna_core.values,immune_core.values]
fused_P, convergence_diffs = snf(p_list, S_list)
#查看convergence_diffs的类型
print(type(convergence_diffs))
取出list里面的第一个元素
print(convergence_diffs[0])


#保存fuse_P
np.savetxt('fused_P.txt', fused_P, fmt='%f')

import matplotlib.pyplot as plt
# Plotting convergence for each data type
plt.figure(figsize=(10, 6))
for i, diffs in enumerate(convergence_diffs):
    plt.plot(diffs, label=f'Data Type {i + 1}')

plt.plot(convergence_diffs[0], label='image data')
plt.plot(convergence_diffs[1], label='mrna data')
plt.plot(convergence_diffs[2], label='mirna data')
plt.plot(convergence_diffs[3], label='immune data')

plt.xlabel('Iteration')
plt.ylabel('Frobenius Norm of Change')
plt.title('Convergence over Iterations')
plt.legend()
plt.show()
plt.savefig('convergence.png',dpi=300)

#对fused_P进行谱聚类，用肘部法找到最佳K
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from sklearn.cluster import SpectralClustering


import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from sklearn.cluster import SpectralClustering

def plot_silhouette_scores(data, max_clusters):
    # 强制矩阵对称
    data = (data + data.T) / 2
    np.fill_diagonal(data, 0)  # 确保对角线为0

    silhouette_scores = []
    for n_clusters in range(2, max_clusters + 1):
        model = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', assign_labels="discretize", random_state=42)
        labels = model.fit_predict(data)
        score = silhouette_score(data, labels, metric='precomputed')
        silhouette_scores.append(score)
        print(f"K={n_clusters} Silhouette Score: {score}")

    plt.figure(figsize=(10, 5))
    plt.plot(range(2, max_clusters + 1), silhouette_scores, marker='o')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Scores by Number of Clusters')
    plt.grid(True)
    plt.show()

# 假设 fused_P 是计算好的相似性矩阵
# 调用函数，最大聚类数设为10
plot_silhouette_scores(fused_P, 10)

#谱聚类聚3类
model = SpectralClustering(n_clusters=3, affinity='precomputed', assign_labels="discretize", random_state=42)
labels = model.fit_predict(fused_P)
#labels的取值各有多少
print(np.bincount(labels))
[212 132 171]
#保存labels
np.savetxt('labels.txt', labels, fmt='%d')

# 检查数据形状
print(mrna.shape)
print(labels.shape)