import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy as np
from matplotlib.colors import ColorConverter
from scipy.stats import ttest_1samp


def getsamplemum(df,sample2index,keyname):
    indexs = df.index
    num_list=[0 for i in range(len(sample2index))]
    for index in indexs:
        num_list[sample2index[index]]=df.loc[index,keyname]
    return num_list

def isbetter(better):
    # 设定假设的比较值，这里假设为 0
    null_hypothesis_value = 0

    # 执行 t 检验
    t_statistic, p_value = ttest_1samp(better, null_hypothesis_value)
    # 打印结果
    print(f"T-statistic: {t_statistic}")
    print(f"P-value: {p_value}")

    # 根据 P-value 判断是否拒绝零假设
    alpha = 0.05
    if p_value < alpha:
        print("拒绝零假设：均值显著大于 0")
    else:
        print("接受零假设：均值不显著大于 0")

def generate_transparent_color(base_color, alpha_factor=0.9):
    # 将颜色从字符串转换为RGB元组
    rgb_base = ColorConverter().to_rgba(base_color)[:3]

    # 生成透明度为原始颜色90%的颜色
    transparent_color = (*rgb_base, alpha_factor)

    return transparent_color



if __name__=="__main__":
    fig=plt.figure(figsize=(1.8, 1.8), dpi=300)
    df=pd.read_csv(r"./realdata/hmBinsPerSample\lorbin.csv",index_col=0)
    index = df.index
    sample2index={}
    for i,sample in enumerate(index):
        sample2index[sample]=i
    basedir = "./realdata/hmBinsPerSample/"
    lorbin = getsamplemum(pd.read_csv(f"{basedir}/lorbin.csv", index_col=0), sample2index, 'm')
    SemiBin = getsamplemum(pd.read_csv(f"{basedir}/semibin.csv", index_col=0), sample2index, 'm')
    SemiBin2 = getsamplemum(pd.read_csv(f"{basedir}/semibin2.csv", index_col=0), sample2index, 'm')
    metabat2 = getsamplemum(pd.read_csv(f"{basedir}/metabat2.csv", index_col=0), sample2index, 'm')
    vamb = getsamplemum(pd.read_csv(f"{basedir}/vamb.csv", index_col=0), sample2index, 'm')
    aamb = getsamplemum(pd.read_csv(f"{basedir}/aamb.csv", index_col=0), sample2index, 'm')
    comebin = getsamplemum(pd.read_csv(f"{basedir}/comebin.csv", index_col=0), sample2index,
                           'm')
    SemiBin_better = sum([1 if (lorbin[i] - SemiBin[i]) >0 else 0 for i in range(len(lorbin))])
    SemiBin2_better = sum([1 if (lorbin[i] - SemiBin2[i]) > 0 else 0 for i in range(len(lorbin))])
    vamb_better = sum([1 if (lorbin[i] - vamb[i]) > 0 else 0 for i in range(len(lorbin))])
    aamb_better = sum([1 if (lorbin[i] - aamb[i]) > 0 else 0 for i in range(len(lorbin))])
    metabat2_better = sum([1 if (lorbin[i] - metabat2[i]) > 0 else 0 for i in range(len(lorbin))])
    comebin_better = sum([1 if (lorbin[i] - comebin[i]) > 0 else 0 for i in range(len(lorbin))])

    SemiBin_worse = sum([-1 if (lorbin[i] - SemiBin[i]) < 0 else 0 for i in range(len(lorbin))])
    SemiBin2_worse = sum([-1 if (lorbin[i] - SemiBin2[i]) < 0 else 0 for i in range(len(lorbin))])
    vamb_worse = sum([-1 if (lorbin[i] - vamb[i]) < 0 else 0 for i in range(len(lorbin))])
    aamb_worse =sum([-1 if (lorbin[i] - aamb[i]) < 0 else 0 for i in range(len(lorbin))])
    metabat2_worse = sum([-1 if (lorbin[i] - metabat2[i]) < 0 else 0 for i in range(len(lorbin))])
    comebin_worse = sum([-1 if (lorbin[i] - comebin[i]) < 0 else 0 for i in range(len(lorbin))])

    SemiBin_same = sum([1 if (lorbin[i] - SemiBin[i]) == 0 else 0 for i in range(len(lorbin))])
    SemiBin2_same = sum([1 if (lorbin[i] - SemiBin2[i]) == 0 else 0 for i in range(len(lorbin))])
    vamb_same = sum([1 if (lorbin[i] - vamb[i]) == 0 else 0 for i in range(len(lorbin))])
    aamb_same = sum([1 if (lorbin[i] - aamb[i]) == 0 else 0 for i in range(len(lorbin))])
    metabat2_same = sum([1 if (lorbin[i] - metabat2[i]) == 0 else 0 for i in range(len(lorbin))])


    x = ['SemiBin', 'SemiBin2', 'VAMB', 'AVAMB', 'MetaBAT2', 'ComeBin']

    data1 = [SemiBin_better, SemiBin2_better, vamb_better, aamb_better, metabat2_better,comebin_better]
    data2 = [SemiBin_worse, SemiBin2_worse, vamb_worse, aamb_worse, metabat2_worse,comebin_worse]

    data3 = [SemiBin_same, SemiBin2_same, vamb_same, aamb_same, metabat2_same]

    print(data1)
    print(data2)

    print('better', [(i / 104.0) * 100 for i in data1])
    print('worse', [(i / 104.0) * 100 for i in data2])
    ax = fig.add_subplot(111)

    import seaborn as sns

    discrete_colors = sns.color_palette('Dark2', n_colors=2)
    worse_color = discrete_colors[0]
    better_color = discrete_colors[1]

    # 绘制好于的水平柱形（向右）
    same_color = generate_transparent_color('steelblue',0.8)
    plt.barh(x, data1, color=worse_color, label='w')
    # 绘制差于的水平柱形（向左）
    plt.barh(x, data2, color=better_color, label='b', left=np.zeros(len(data2)))

    for data in data1:
        print(data)


    for data in data2:
        print(data)

    # 设置x轴刻度
    # plt.xticks(np.arange(min(min(data1), min(data2)), max(max(data1), max(data2)) + 1, 2))
    #plt.xticks([-15, 0, 20, 50])

    # 添加0刻度线
    plt.axvline(0, color='black', linewidth=0.8)

    # 添加标签和标题

    # 添加图例
    plt.legend(ncol=1,frameon=False,fontsize=7)
    plt.xticks([-30,0,30,90],['30','0','30','90'],fontsize=8)
    plt.yticks('')
    plt.xlabel('# sample',fontsize=8)
    plt.title('mBins',fontsize=8)
    ax.patch.set_facecolor('none')
    fig.set_facecolor('none')

    # 显示图形
    plt.savefig(f'{basedir}/mMAG_%_.svg', bbox_inches='tight', pad_inches=0.02,format='svg')
    plt.show()

    """

    isbetter(SemiBin)
    isbetter(SemiBin2)
    isbetter(vamb)
    isbetter(aamb)
    isbetter(metabat2)



    data = [SemiBin,SemiBin2,vamb,aamb,metabat2]
    boxprops = dict(linewidth=1)
    flierprops = dict(marker='o', markersize=3)

    # 绘制箱线图
    plt.boxplot(data,vert=False,labels=['SemiBin','SemiBin2','VAMB','AVAMB','MetaBAT2',],boxprops=boxprops,flierprops=flierprops)

    # 添加标题和标签
    plt.title('')
    plt.xlabel('# mMAGs LorBin better per sample',fontsize=8)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)

    # 显示图形

    plt.savefig('#hmMAGpersample_better_m_.png', bbox_inches='tight', pad_inches=0.02)
    plt.show()
    """
