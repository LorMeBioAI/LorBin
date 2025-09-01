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
base_color = "#008080"

# 生成更深和更浅的颜色

if __name__=="__main__":
    #x=['airways','GI','oral','skin','urog']

    color1 = generate_transparent_color(base_color, 0.8)  # 调整亮度为 0.7，可以根据需要调整
    color2 = generate_transparent_color(base_color, 0.6)  # 调整亮度为 1.3，可以根据需要调整
    color3 = generate_transparent_color(base_color, 0.4)  # 调整亮度为 1.5，可以根据需要调整
    color4 = generate_transparent_color(base_color, 0.3)  # 调整亮度为 1.5，可以根据需要调整
    color5 = generate_transparent_color(base_color, 0.2)  # 调整亮度为 1.5，可以根据需要调整
    #color6 = generate_transparent_color(base_color, 0.4)  # 调整亮度为 1.5，可以根据需要调整

    airways = pd.read_csv(r"./CAMI/bins/airways.csv",index_col=0)
    #x = ['LorBin','SemiBin2','AVAMB','SemiBin2','VAMB','MetaBAT2']
    x=['LorBin','COMEBin','SemiBin','SemiBin2','VAMB','AAMB','Metabat2']
    """
    LorBin=[139, 171, 226,186,126]
    VAMB=[63,86,124,72,78]
    AAMB=[74,98,118,90,70]
    SemiBin=[96,140,159,138,112]
    SemiBin2=[94,149,170,141,118]
    """
    #Metabat2=[90,120,136,113,104]
    airways_h=list(airways.loc[:,'hMAGs'])
    airways_90 = list(airways.loc[:,'MAGs90'])
    airways_80 = list(airways.loc[:,'MAGs80'])
    airways_70 = list(airways.loc[:,'MAGs70'])
    airways_60 = list(airways.loc[:,'MAGs60'])
    airways_50 =list(airways.loc[:,'MAGs50'])

    gi = pd.read_csv('./CAMI/bins//gi.csv',index_col=0)
    gi_h = list(gi.loc[:,'hMAGs'])
    gi_90 = list(gi.loc[:,'MAGs90'])
    gi_80 = list(gi.loc[:,'MAGs80'])
    gi_70 = list(gi.loc[:,'MAGs70'])
    gi_60 = list(gi.loc[:,'MAGs60'])
    gi_50 = list(gi.loc[:,'MAGs50'])

    oral = pd.read_csv("./CAMI/bins//oral.csv",index_col=0)
    oral_h = list(oral.loc[:,'hMAGs'])
    oral_90 = list(oral.loc[:,'MAGs90'])
    oral_80 = list(oral.loc[:,'MAGs80'])
    oral_70 = list(oral.loc[:,'MAGs70'])
    oral_60 = list(oral.loc[:,'MAGs60'])
    oral_50 = list(oral.loc[:,'MAGs50'])

    skin = pd.read_csv("./CAMI/bins//skin.csv",index_col=0)
    skin_h = list(skin.loc[:,'hMAGs'])
    skin_90 = list(skin.loc[:,'MAGs90'])
    skin_80 = list(skin.loc[:,'MAGs80'])
    skin_70 = list(skin.loc[:,'MAGs70'])
    skin_60 = list(skin.loc[:,'MAGs60'])
    skin_50 = list(skin.loc[:,'MAGs50'])

    urog = pd.read_csv('./CAMI/bins/urog.csv',index_col=0)
    urog_h = list(urog.loc[:,'hMAGs'])
    urog_90 = list(urog.loc[:,'MAGs90'])
    urog_80 = list(urog.loc[:,'MAGs80'])
    urog_70 = list(urog.loc[:,'MAGs70'])
    urog_60 = list(urog.loc[:,'MAGs60'])
    urog_50 = list(urog.loc[:,'MAGs50'])

    fig=plt.figure(figsize=(6.8,1.2),dpi=300)



    ax=plt.subplot(1,5,1)
    plt.title('airways',fontsize=8)
    plt.barh(x, airways_h, height=0.6, label='>90%,cont<5%', color=base_color)
    plt.barh(x, airways_90, height=0.6,left=airways_h,label='>90%',color=color1)
    plt.barh(x, airways_80, height=0.6,left=[airways_h[i]+airways_90[i] for i in range(7)], label='>80%', color=color2)
    plt.barh(x, airways_70, height=0.6,left=[airways_h[i]+airways_90[i]+airways_80[i] for i in range(7)], label='>70%', color=color3)
    plt.barh(x, airways_60, height=0.6,left=[airways_h[i]+airways_90[i]+airways_80[i]+airways_70[i] for i in range(7)], label='>60%', color=color4)
    plt.barh(x, airways_50, height=0.6,
             left=[airways_h[i] + airways_90[i] + airways_80[i] + airways_70[i] + airways_60[i] for i in range(7)], label='>50%',
             color=color5)
    ax.patch.set_facecolor('none')
    ax.legend(loc='lower center', bbox_to_anchor=(2.75, -0.45), fancybox=True, shadow=True, ncol=7,
               frameon=False, fontsize=7)

    """
    plt.barh(x, airways_40, height=0.6,
             left=[airways_h[i] + airways_90[i] + airways_80[i] + airways_70[i]+airways_60[i]+airways_50[i] for i in range(6)], label='>40%',
             color=color6)
    """
    plt.yticks(x,fontsize=8)
    plt.xticks(fontsize=8)

    ax=plt.subplot(1, 5, 2)
    plt.title('gi',fontsize=8)
    plt.barh(x, gi_h, height=0.6,label='com>90%,con<5%',color=base_color)
    plt.barh(x, gi_90, left = gi_h, height=0.6,label='>90%',color=color1)
    plt.barh(x, gi_80, height=0.6, left=[gi_h[i]+gi_90[i] for i in range(7)], label='>80%', color=color2)
    plt.barh(x, gi_70, height=0.6, left=[gi_h[i]+gi_90[i] + gi_80[i] for i in range(7)], label='>70%',
             color=color3)
    plt.barh(x, gi_60, height=0.6, left=[gi_h[i]+gi_90[i] + gi_80[i] + gi_70[i] for i in range(7)],
             label='>60%', color=color4)
    plt.barh(x, gi_50, height=0.6, left=[gi_h[i]+gi_90[i] + gi_80[i] + gi_70[i]+gi_60[i] for i in range(7)],
             label='>50%', color=color5)
    #plt.barh(x, gi_40, height=0.6, left=[gi_h[i] + gi_90[i] + gi_80[i] + gi_70[i]+gi_60[i]+gi_50[i] for i in range(6)],label='>40%', color=color6)
    plt.yticks([])
    plt.xticks(fontsize=8)
    ax.patch.set_facecolor('none')

    ax=plt.subplot(1, 5, 3)
    plt.title('oral',fontsize=8)
    plt.barh(x, oral_h, height=0.6,label='com>90%,con<5%',color=base_color)
    plt.barh(x, oral_90,left=oral_h, height=0.6,label='>90%',color=color1)
    plt.barh(x, oral_80, height=0.6, left=[oral_h[i]+oral_90[i] for i in range(7)], label='>80%', color=color2)
    plt.barh(x, oral_70, height=0.6, left=[oral_h[i]+oral_90[i] + oral_80[i] for i in range(7)], label='>70%',
             color=color3)
    plt.barh(x, oral_60, height=0.6, left=[oral_h[i]+oral_90[i] + oral_80[i] + oral_70[i] for i in range(7)],
             label='>60%', color=color4)
    plt.barh(x, oral_50, height=0.6, left=[oral_h[i] + oral_90[i] + oral_80[i] + oral_70[i] +oral_60[i]for i in range(7)],
             label='>50%', color=color5)
    #plt.barh(x, oral_40, height=0.6, left=[oral_h[i] + oral_90[i] + oral_80[i] + oral_70[i] +oral_60[i]+oral_50[i]for i in range(6)],label='>40%', color=color6)
    plt.yticks([])
    plt.xticks(fontsize=8)
    ax.patch.set_facecolor('none')

    ax=plt.subplot(1, 5, 4)
    plt.title('skin',fontsize=8)
    plt.barh(x, skin_h, height=0.6,label='com>90%,con<5%',color=base_color)
    plt.barh(x, skin_90, left=skin_h, height=0.6,label='>90%',color=color1)
    plt.barh(x, skin_80, height=0.6, left=[skin_h[i]+skin_90[i] for i in range(7)], label='>80%', color=color2)
    plt.barh(x, skin_70, height=0.6, left=[skin_h[i]+skin_90[i] + skin_80[i] for i in range(7)], label='>70%',
             color=color3)
    plt.barh(x, skin_60, height=0.6, left=[skin_h[i]+skin_90[i] + skin_80[i] + skin_70[i] for i in range(7)],
             label='>60%', color=color4)
    plt.barh(x, skin_50, height=0.6, left=[skin_h[i]+skin_90[i] + skin_80[i] + skin_70[i] +skin_60[i] for i in range(7)],
             label='>50%', color=color5)
    #plt.barh(x, skin_40, height=0.6, left=[skin_h[i]+skin_90[i] + skin_80[i] + skin_70[i] +skin_60[i]+skin_50[i] for i in range(6)],label='>40%', color=color6)
    plt.yticks([])
    plt.xticks(fontsize=8)
    ax.patch.set_facecolor('none')

    ax=plt.subplot(1, 5, 5)
    plt.title('urog',fontsize=8)
    plt.barh(x, urog_h, height=0.6,label='com>90%,con<5%',color=base_color)
    plt.barh(x, urog_90,left=urog_h, height=0.6,label='>90%',color=color1)
    plt.barh(x, urog_80, height=0.6, left=[urog_h[i]+urog_90[i]for i in range(7)], label='>80%', color=color2)
    plt.barh(x, urog_70, height=0.6, left=[urog_h[i]+urog_90[i] + urog_80[i] for i in range(7)], label='>70%',
             color=color2)
    plt.barh(x, urog_60, height=0.6, left=[urog_h[i]+urog_90[i] + urog_80[i] + urog_70[i] for i in range(7)],
             label='>60%', color=color4)
    plt.barh(x, urog_50, height=0.6, left=[urog_h[i] + urog_90[i] + urog_80[i] + urog_70[i]+urog_60[i] for i in range(7)],
             label='>50%', color=color5)
    #plt.barh(x, urog_40, height=0.6, left=[urog_h[i] + urog_90[i] + urog_80[i] + urog_70[i]+urog_60[i]+urog_50[i] for i in range(6)],label='>40%', color=color6)
    plt.yticks([])
    plt.xticks(fontsize=8)
    ax.patch.set_facecolor('none')

    fig.set_facecolor('none')

    plt.figtext(0.35,-0.3,'# Bins(< 10% contamination) ',fontsize=8)

    # 调整子图的布局
    #plt.subplots_adjust(bottom=-0.3)

    # 显示图形
    plt.savefig('CAMI_MAGs.svg', bbox_inches='tight', pad_inches=0.02,format='svg')
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