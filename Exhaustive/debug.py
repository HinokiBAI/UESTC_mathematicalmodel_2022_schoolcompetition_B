# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import Normalizer
import matplotlib

# %%
person = 'A'
reference = pd.read_csv('./data/Refrence.csv')
# 计算每个个体的基因座上的等位基因
A_refrence_data = reference[reference['sample'] == person]
alleles_dict = {A_refrence_data['marker']:[A_refrence_data['allele1'],A_refrence_data['allele2']] for index,A_refrence_data in A_refrence_data.iterrows()}
delete_list = ['Yindel', 'AMEL', 'DYS391']

# %%
init_data = pd.read_csv(f'data/{person}100pg.hid_Genotype.csv', usecols=['Marker', 'Allele', 'Height', 'Size'])
init_data = init_data[(init_data['Size'] != 'Dropout') & ~(init_data['Marker'].isin(delete_list))]


# %%
threshold = 5200
# 阈值列表
threshold_list = []
# 剩余噪声数目
remain_noise_list = []
# 丢失基因数目
lost_gene_list = []

test = None
while threshold <= 5200:
    error = 0
    remain_gene = 0
    # 过滤掉的数据
    filter_data = init_data[init_data['Height'] < threshold]
    # 过滤后的数据
    remain_data = init_data[init_data['Height'] >= threshold]
    threshold_list.append(threshold)

    for index, row in filter_data.iterrows():
        if row['Allele'] == str(alleles_dict[row['Marker']][0]) or row['Allele'] == str(alleles_dict[row['Marker']][1]):
            error += 1

    lost_gene_list.append(error)

    for index, row in remain_data.iterrows():
        if row['Allele'] != alleles_dict[row['Marker']][0] and row['Allele'] != alleles_dict[row['Marker']][1]:
            remain_gene += 1
    remain_noise_list.append(remain_gene)
    threshold += 0.5


# %%
def f(x, y):
    return np.sqrt(x ** 2 + y ** 2) *100

# %%
x_norm = Normalizer(norm='max').fit_transform(np.array(remain_noise_list, dtype='float32').reshape(1, -1))
y_norm = Normalizer(norm='max').fit_transform(np.array(lost_gene_list, dtype='float32').reshape(1, -1))

area_list = [f(x, y) for x, y in zip(x_norm, y_norm)]
matplotlib.rcParams['font.sans-serif'] = ['SimHei']
matplotlib.rcParams['font.family']='sans-serif' #解决负号'-'显示为方块的问题
matplotlib.rcParams['axes.unicode_minus'] = False
plt.figure(dpi=300)
plt.title(f'{person}基因剩余噪声数目与丢失基因数目的关系')
plt.xlabel('剩余噪声数目')
plt.ylabel('丢失基因数目')
plt.scatter(x_norm, y_norm, marker='o', alpha=0.2, s=area_list, c='c')

min_area = min(area_list[0])
min_index = area_list[0].tolist().index(min_area)
# print(min_area, min_index)
min_x = x_norm[0][min_index]
min_y = y_norm[0][min_index]
plt.scatter(min_x, min_y, marker='o', alpha=1, c='r')
print(f'阈值为{threshold_list[min_index]}')
print(f'剩余噪声数目为{remain_noise_list[min_index]}')
print(f'丢失基因数目为{lost_gene_list[min_index]}')


