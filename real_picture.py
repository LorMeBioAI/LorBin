import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb, to_hex, ColorConverter


def generate_transparent_color(base_color, alpha_factor=0.9):
    # 将颜色从字符串转换为RGB元组
    rgb_base = ColorConverter().to_rgba(base_color)[:3]

    # 生成透明度为原始颜色90%的颜色
    transparent_color = (*rgb_base, alpha_factor)

    return transparent_color

# 定义基准颜色（steelblue）
base_color = "#4682B4"
color = generate_transparent_color(base_color, 0.8)
# 创建示例数据
categories = ['LorBin', 'SemiBin','SemiBin2','COMEBin','VAMB','AAMB','MetaBAT2']
x_len = np.arange(len(categories))
total_width, n = 0.6, 1
width = 0.6
fig=plt.figure(figsize=(1.8, 1.8), dpi=300)
values1 = [224,129,202,140,74,107,80]
values2 = [455,414,469,327,209,284,80]
ax =fig.add_subplot(111)

import seaborn as sns

discrete_colors = sns.color_palette('summer', n_colors=2)
base_color=discrete_colors[0]
color=discrete_colors[1]

# 画柱状堆积图
plt.barh(categories,values1, width,align='center', color=base_color,tick_label=categories,label='h')
plt.barh(categories, values2, width,left=values1, label='m',color=color,tick_label=categories)
plt.yticks(categories,fontsize=8)
plt.xticks(fontsize=8)
yticks = x_len - width

plt.text(values1[0]+values2[0]+40, yticks[0]+0.5, values1[0]+values2[0], ha='center',  fontsize=6,
                 zorder=10)
plt.text(values1[1]+values2[1]+40, yticks[1]+0.5, values1[1]+values2[1], ha='center',  fontsize=6,
                 zorder=10)
plt.text(values1[2]+values2[2]+40, yticks[2]+0.5, values1[2]+values2[2], ha='center',  fontsize=6,
                 zorder=10)
plt.text(values1[3]+values2[3]+40, yticks[3]+0.5, values1[3]+values2[3], ha='center',  fontsize=6,
                 zorder=10)
plt.text(values1[4]+values2[4]+40, yticks[4]+0.5, values1[4]+values2[4], ha='center',  fontsize=6,
                 zorder=10)
plt.text(values1[5]+values2[5]+40, yticks[5]+0.5, values1[5]+values2[5], ha='center',  fontsize=6,
                 zorder=10)
plt.text(values1[6]+values2[6]+40, yticks[6]+0.5, values1[6]+values2[6], ha='center',  fontsize=6,
                 zorder=10)



#plt.ylabel('# of bins')
plt.xlim(0,720)
#plt.xticks([0,250,500,800])
plt.xticks([0,300,600,900])
plt.xlabel('# hmBins',fontsize=8)

# 添加图例
plt.legend(prop={'size':6},ncol=1,frameon=False)

# 显示图形
fig.set_facecolor('none')
ax.patch.set_facecolor('none')
plt.savefig('realdataset_hifiasm.svg', bbox_inches='tight', pad_inches=0.02,format='svg')
plt.show()
