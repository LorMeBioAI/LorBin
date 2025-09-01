import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib.colors import to_hex, hex2color, rgb_to_hsv, hsv_to_rgb, to_rgb, ColorConverter


def adjust_brightness(color, factor):
    rgb_color = np.clip(np.array(to_rgb(color)) * factor, 0, 1)
    return to_hex(rgb_color)

def generate_transparent_color(base_color, alpha_factor=0.9):
    # 将颜色从字符串转换为RGB元组
    rgb_base = ColorConverter().to_rgba(base_color)[:3]

    # 生成透明度为原始颜色90%的颜色
    transparent_color = (*rgb_base, alpha_factor)

    return transparent_color

# 定义基准颜色（steelblue）
base_color = "#4682B4"

# 生成更深和更浅的颜色

if __name__=="__main__":
    #x=['airways','GI','oral','skin','urog']

    other_color=generate_transparent_color(base_color,0.75)

    x=['LorBin','COMEBin','SemiBin','SemiBin2','VAMB','AAMB','Metabat2']
    """
    LorBin=[139, 171, 226,186,126]
    VAMB=[63,86,124,72,78]
    AAMB=[74,98,118,90,70]
    SemiBin=[96,140,159,138,112]
    SemiBin2=[94,149,170,141,118]
    """
    #Metabat2=[90,120,136,113,104]
    df = pd.read_excel('./CAMI/evaluation/acc.xlsx',index_col=0)
    airways=list(df.loc['airways',:])
    gi = list(df.loc['gi',:])
    oral = list(df.loc['oral',:])
    skin =list(df.loc['skin',:])
    urog = list(df.loc['urog',:])


    fig=plt.figure(figsize=(6.8,1.2),dpi=300)



    SemiBin_color= "#b0b1b6"
    SemiBin2_color = "#e4e6e1"
    AAMB_color ="#c2cedc"
    LorBin_color ="#d4baad"
    VAMB_color ="#e1ccb1"

    discrete_colors = sns.color_palette('Set2', n_colors=len(x))

    ax=plt.subplot(1,5,1)
    #plt.title('airways',fontsize=8)
    plt.barh(x, airways, height=0.6,color=discrete_colors)
    plt.yticks(x,fontsize=8)
    plt.xticks(fontsize=8)
    ax.patch.set_facecolor('none')

    ax=plt.subplot(1, 5, 2)
    #plt.title('gi',fontsize=8)
    plt.barh(x, gi, height=0.6,color=discrete_colors)
    plt.yticks([])
    plt.xticks(fontsize=8)
    ax.patch.set_facecolor('none')

    ax=plt.subplot(1, 5, 3)
    #plt.title('oral',fontsize=8)
    plt.barh(x, oral, height=0.6,color=discrete_colors)
    plt.yticks([])
    plt.xticks(fontsize=8)
    ax.patch.set_facecolor('none')

    ax=plt.subplot(1, 5, 4)
    #plt.title('skin',fontsize=8)
    plt.barh(x, skin, height=0.6,color=discrete_colors)
    plt.yticks([])
    plt.xticks(fontsize=8)
    ax.patch.set_facecolor('none')

    ax=plt.subplot(1, 5, 5)
    #plt.title('urog',fontsize=8)
    plt.barh(x, urog, height=0.6,color=discrete_colors)
    plt.yticks([])
    plt.xticks(fontsize=8)
    plt.figtext(0.45,-0.2,'Accuracy ',fontsize=8)
    ax.patch.set_facecolor('none')
    fig.set_facecolor('none')
    """
    fig.legend(loc='lower center', bbox_to_anchor=(1.6,-0.45), fancybox=True, shadow=True, ncol=4,
               frameon=False,fontsize=7)
    """

    # 调整子图的布局
    #plt.subplots_adjust(bottom=-0.3)

    # 显示图形
    plt.savefig('CAMI_MAGs_acc.svg', bbox_inches='tight', pad_inches=0.02,format='svg')
    plt.show()

    """

    bar_width = 0.2
    index = np.arange(len(x))

    # 创建横向柱状图
    plt.barh(index, LorBin, height=bar_width, label='LorBin')
    plt.barh(index, VAMB, height=bar_width, label='VAMB')
    plt.barh(index, AAMB, height=bar_width, label='AAMB', left=[LorBin[i] + VAMB[i] for i in range(len(LorBin))])
    plt.barh(index, SemiBin, height=bar_width, label='SemiBin',
             left=[LorBin[i] + VAMB[i] + AAMB[i] for i in range(len(LorBin))])
    plt.barh(index, SemiBin2, height=bar_width, label='SemiBin2',
             left=[LorBin[i] + VAMB[i] + AAMB[i] + SemiBin[i] for i in range(len(LorBin))])

    # 添加标签和标题
    plt.xlabel('Values')
    plt.ylabel('Categories')
    plt.title('Horizontal Bar Chart for Categories')
    plt.yticks(index, x)
    plt.legend()

    # 显示图表
    plt.show()
    """