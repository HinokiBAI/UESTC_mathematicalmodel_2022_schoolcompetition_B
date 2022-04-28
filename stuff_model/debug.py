# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from scipy.optimize import curve_fit
import scipy.optimize as optimize
from sklearn.linear_model import LinearRegression

# %%
persons = ['A', 'B', 'C']
nums = [10, 25, 50, 100, 200, 500]
# nums = [10]

delete_list = ['Yindel', 'AMEL', 'DYS391', 'Amel']
reference = pd.read_csv('data/Refrence.csv')
# 计算每个个体的基因座上的等位基因
alleles_dict_list = []
for person in persons:
    A_refrence_data = reference[reference['sample'] == person]
    A_refrence_data = A_refrence_data[~A_refrence_data['marker'].isin(delete_list)]
    alleles_dict = {A_refrence_data['marker']:[A_refrence_data['allele1'],A_refrence_data['allele2']] for index,A_refrence_data in A_refrence_data.iterrows()}
    alleles_dict_list.append(alleles_dict)


# %%
goal_data = {'基因座':[], '斜率':[], '截距':[], 'R2':[]}
radio_gate = 0.5


marker = list(alleles_dict_list[0].keys())[1] # *

for marker in alleles_dict_list[0].keys():
    matplotlib.rcParams['font.sans-serif'] = ['SimHei']
    matplotlib.rcParams['font.family']='sans-serif' #解决负号'-'显示为方块的问题
    matplotlib.rcParams['axes.unicode_minus'] = False
    plt.figure(dpi=300)
    plt.title(f'{marker}基因座的Allele与隐峰比例拟合指数曲线')
    plt.xlabel('Allele')
    plt.ylabel('隐峰比例')
    # plt.ylim(-0.5, 1)
    all_data_dict = {}

    for index, person in enumerate(['A', 'B', 'C']):
        radio_list_list = []
        x_list = []
        y_list = []
        for num in nums:
            init_data = pd.read_csv(f'data/单个个体/{person}{num}pg.hid_Genotype.csv', usecols=['Marker', 'Allele', 'Height', 'Size'])
            init_data = init_data[(init_data['Size'] != 'Dropout') & ~(init_data['Marker'].isin(delete_list))]
            init_data['Height'] = init_data['Height'].apply(pd.to_numeric, errors='coerce').fillna(0.0)
            init_data['Allele'] = init_data['Allele'].apply(pd.to_numeric, errors='coerce').fillna(0.0)
            # 获得每个基因座的数据
            marker_data = init_data[init_data['Marker'] == marker]
            # 获得每个基因座的所有等位基因
            alleles_list = marker_data['Allele'].values
            # 获得真实等位基因
            true_alleles = np.array(alleles_dict_list[index][marker], dtype=float) # *

            # 获得等位基因横坐标
            x = []
            # 获得真实高度
            if true_alleles[0] == true_alleles[1]:
                for i in alleles_list:
                    if i < true_alleles[0]:
                        x.append(i)
                true_height = marker_data[marker_data['Allele'] == true_alleles[0]]['Height'].values[0] / 2
                height_list = marker_data[marker_data['Allele'].isin(x)]['Height'].values
                radio_list = height_list / true_height

            if true_alleles[0] != true_alleles[1]:
                for i in alleles_list:
                    if i < true_alleles[1] and  i != true_alleles[0]:
                        x.append(i)
                true_height = (marker_data[marker_data['Allele'] == true_alleles[0]]['Height'].values[0] + marker_data[marker_data['Allele'] == true_alleles[1]]['Height'].values[0]) / 2
                height_list = marker_data[marker_data['Allele'].isin(x)]['Height'].values
                radio_list = height_list / true_height

            flag = []
            for i in range(len(radio_list)):
                if radio_list[i] > radio_gate:
                    flag.append(i)
            radio_list = np.delete(radio_list, flag)
            radio_list_list.append(radio_list)
            x = np.delete(x, flag)
            plt.scatter(x, radio_list, color='blue', s=12)
            for i in range(len(x)):
                # 判断是否有all_data_dict[x[i]]，若有，则追加，若没有则创建
                if x[i] in all_data_dict.keys():
                    all_data_dict[x[i]].append(radio_list[i])
                else:
                    all_data_dict[x[i]] = [radio_list[i]]


    # 绘制all_data_dict中的数据
    for key in all_data_dict.keys():
        x_list.append(key)
        y_list.append(np.mean(all_data_dict[key]))
    plt.scatter(x_list, y_list, color='red', s=12)
    model = LinearRegression()
    x = np.array(x_list).reshape(-1, 1)
    y = np.array(y_list).reshape(-1, 1)
    model.fit(x, y)
    plt.plot(x, model.predict(x), color='#dd3032', linewidth=1)
    goal_data['基因座'].append(marker)
    goal_data['斜率'].append(model.coef_[0][0])
    goal_data['截距'].append(model.intercept_[0])
    goal_data['R2'].append(model.score(x, y))




