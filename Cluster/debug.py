# %%
import pandas as pd

# %%
# 数据初始化
person = 'A'
num = 10
category = '单个个体'

# 阈值
threshold = 25

delete_list = ['Yindel', 'AMEL', 'DYS391', 'Amel']
reference = pd.read_csv('data/Refrence.csv')
# 计算每个个体的基因座上的等位基因
A_refrence_data = reference[reference['sample'] == person]
A_refrence_data = A_refrence_data[~A_refrence_data['marker'].isin(delete_list)]
alleles_dict = {A_refrence_data['marker']:[A_refrence_data['allele1'],A_refrence_data['allele2']] for index,A_refrence_data in A_refrence_data.iterrows()}


init_data = pd.read_csv(f'data/{category}/{person}{num}pg.hid_Genotype.csv', usecols=['Marker', 'Allele', 'Height', 'Size'])
init_data = init_data[(init_data['Size'] != 'Dropout') & ~(init_data['Marker'].isin(delete_list))]
init_data['Height'] = init_data['Height'].apply(pd.to_numeric, errors='coerce').fillna(0.0)

stutter = {}
# 遍历alleles_dict
for key, value in alleles_dict.items():
    new_list = []
    for i in value:
        new_list.append(float(i)-1)
    stutter[key] = new_list

# 阴性对照
negative_data = pd.read_csv('data/阴性对照NC.hid_Genotype.csv', usecols=['Marker', 'Allele', 'Height', 'Size'])
negative_data = negative_data[(negative_data['Size'] != 'Dropout') & ~(negative_data['Marker'].isin(delete_list))]
negative_data['Height'] = negative_data['Height'].apply(pd.to_numeric, errors='coerce').fillna(0.0)

negative_data.head()

# %%
# 1. 噪声过滤
remain_data = init_data[init_data['Height'] >= threshold]

# 2. 真实基因剔除
for key, value in alleles_dict.items():
    remain_data = remain_data.drop(remain_data[(remain_data['Marker'] == key) & ((remain_data['Allele'] == float(value[0])) | (remain_data['Allele'] == float(value[1])))].index)

print(remain_data)
print('--------------------')
# 3. 删去真实等位基因 -1 的峰
for key, value in stutter.items():
    remain_data = remain_data.drop(remain_data[(remain_data['Marker'] == key) & ((remain_data['Allele'] == float(value[0])) | (remain_data['Allele'] == float(value[1])))].index)

print(remain_data)
print('--------------------')

# 4. 减去阴性对照峰值
# 遍历remain_data
for index, row in remain_data.iterrows():
    for key, value in negative_data.iterrows():
        if row['Marker'] == value['Marker']:
            if row['Allele'] == value['Allele']:
                remain_data.loc[index, 'Height'] = row['Height'] - value['Height']




