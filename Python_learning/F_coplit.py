######################################################################
# xingjinyin 
# 2024-12-22
# 使用常用的python 内部数据库，画一个渐变小提琴图
####################################################################
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
# 生成数据
np.random.seed(0)
data = np.random.normal(loc=0, scale=1, size=100)
data2 = np.random.normal(loc=0, scale=1, size=100)
data3 = np.random.normal(loc=0, scale=1, size=100)
data4 = np.random.normal(loc=0, scale=1, size=100)
data5 = np.random.normal(loc=0, scale=1, size=100)
data6 = np.random.normal(loc=0, scale=1, size=100)

# 画图并保存
plt.figure(figsize=(8, 6))
# 颜色渐变
sns.set_palette('coolwarm')
# 画出小提琴图
sns.violinplot(data=[data, data2, data3, data4], linewidth=2, palette='coolwarm')
# 画出均值
plt.plot([0, 1, 2, 3], [np.mean(data), np.mean(data2), np.mean(data3),np.mean(data4)], 'o', color='black')
# 画出标准差
plt.plot([0, 1, 2, 3], [np.mean(data) + np.std(data), np.mean(data2) + np.std(data2), np.mean(data3) + np.std(data3),np.mean(data4) + np.std(data4)], '--', color='black')
# 把图像背景变成空白
plt.gca().patch.set_facecolor('white')
# 图像背景透明度
plt.gca().patch.set_alpha(0)
# 字体大小
plt.rcParams['font.size'] = 16
plt.title('Violin Plot')
# 横纵坐标图例字体大小
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
# x轴y轴标签大小
plt.xlabel('Data', fontsize=16)
plt.ylabel('Value', fontsize=16)
plt.grid(True)
plt.savefig('violin_plot.png')
plt.show()
























