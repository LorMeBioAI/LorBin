import plotly.graph_objects as go
import pandas as pd
import os
from Bio import SeqIO
from matplotlib.colors import ColorConverter
from matplotlib import colors


def generate_transparent_color(base_color, alpha_factor=0.9):
    # 将颜色从字符串转换为RGB元组
    rgb_base = ColorConverter().to_rgba(base_color)[:3]

    # 生成透明度为原始颜色90%的颜色
    transparent_color = (*rgb_base, alpha_factor)

    return transparent_color

# 定义节点信息
nodes = ["all contigs", "first clustering","noises and poor bins","second clustering","noises",'hmMAGs','lMAGs']
df=pd.read_csv("./realdata/lorbin_sankey.csv",index_col=0)
df1 = df[df['stage']==1]
df2 = df[df['stage']==2]
df1_hm = df[(df['stage']==1)&(df['quality']!='l')]
bin1 = list(df1.loc[:,'bin_size'])
bin1_hm = list(df1_hm.loc[:,'bin_size'])
print(sum(bin1_hm))
bin2 = list(df2.loc[:,'bin_size'])


first_bin = sum(bin1)
total_contigs=403608
recluster = total_contigs - first_bin
noise = 403608 - 278342
second_bin_hm = sum(bin2)
second_bin_l = 278342 - second_bin_hm - first_bin
first_bin_l = first_bin - sum(bin1_hm)

"""
flowcolor1= colors.to_rgba('red')[:3]
flowcolor1 = tuple(int(component*255+0.5) for component in flowcolor1)
flowcolor1 = flowcolor1+(0.2,)

flowcolor2= colors.to_rgba('orange')[:3]
flowcolor2 = tuple(int(component*255+0.5) for component in flowcolor2)
flowcolor2 = flowcolor2+(0.2,)

flowcolor3= colors.to_rgba('blue')[:3]
flowcolor3 = tuple(int(component*255+0.5) for component in flowcolor3)
flowcolor3 = flowcolor3+(0.2,)

flowcolor4= colors.to_rgba('skyblue')[:3]
flowcolor4 = tuple(int(component*255+0.5) for component in flowcolor4)
flowcolor4 = flowcolor4+(0.2,)
"""
flowcolor1 =(150,202,193,0.5)
flowcolor2 =(193,190,214,0.5)
flowcolor3 =(138,175,201,0.5)
flowcolor4 =(193,190,214,0.5)

import seaborn as sns
discrete_colors = sns.color_palette('Dark2', n_colors=7)

links_color =['#96cac1','#c1bed6','#c1bed6','#8aafc9','#8aafc9','#afcf78','#f5f5bd']




links =[{'source':0,'target':1,'value':first_bin},{'source':0,'target':2,'value':recluster},
        {'source':2,'target':3,'value':second_bin_l+second_bin_hm},{'source':2,'target':4,'value':noise},
        {'source':1,'target':5,'value':sum(bin1_hm)},{'source':1,'target':6,'value':first_bin_l},
        {'source':3,'target':5,'value':second_bin_hm},{'source':3,'target':6,'value':second_bin_l}]


# 创建桑基图
fig = go.Figure(go.Sankey(
    node=dict(
        pad=10,
        thickness=15,
        line=dict(color="gray", width=0.15),
        color=links_color,
        #label=nodes,
        #customdata=[link["value"] for link in links]
    ),
    link=dict(
        source=[link["source"] for link in links],
        target=[link["target"] for link in links],
        value=[link["value"] for link in links],
        #color="lightgray"
        #opacity=0.2,
        color=[f'rgba{flowcolor1}',f'rgba{flowcolor1}',f'rgba{flowcolor2}',f'rgba{flowcolor2}',
               f'rgba{flowcolor4}',f'rgba{flowcolor4}',f'rgba{flowcolor3}',f'rgba{flowcolor3}'],

    )
))

# 设置图形布局
fig.update_layout(title="binning",
                  font_size=8, width=300, height=350)

# 显示图形
fig.show()
