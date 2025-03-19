import pandas as pd
import numpy as np
# 数据导入
df = pd.read_csv("Raw_data.csv",header=16) #不轻易改原始文件
# print(df)
# 数据处理和统计
df['金额']=pd.to_numeric(df['金额(元)'].str[1:],errors='ignore')
df['交易时间']=pd.to_datetime(df['交易时间'],errors='ignore')
df['时间']=df['交易时间'].apply(lambda x:x.strftime('%H'))
df['月份']=df['交易时间'].apply(lambda x:x.strftime('%Y-%m'))

#保存文件
# outfile = 'out.csv'
# df.to_csv(outfile)
# print(df)
# 交易时间
df2 = pd.pivot_table(df,index= ['月份'],values = ['金额'],aggfunc= (sum,len,max,np.mean))
print(df2)

#交易对方
df3 = pd.pivot_table(df,index= ['交易对方'],values = ['金额','收/支'],aggfunc= (sum,len))
print(df3)
# 数据可视化







