import statistics

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


if __name__=="__main__":
    envs=['airways','gi','oral','skin','urog']
    fig = plt.figure(figsize=(6.8, 1.2), dpi=300)
    x = ['LorBin','COMEBin','SemiBin','SemiBin2','VAMB','AAMB','MetaBAT2']

    """
    plt.rcParams['font.sans-serif'] = ['SimSun']  # Windows用SimSun，macOS/Linux用'STSong'
    plt.rcParams['font.size'] = 12  # 五号字约为12pt
    plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问
    """

    discrete_colors = sns.color_palette('Set2', n_colors=len(x))
    for n,env in enumerate(envs):
        df = pd.read_csv(f"./unbalance/area_{env}.csv",index_col=0)
        ax = plt.subplot(1,5,n+1)
        y=[]
        y_e=[]
        lorbin = df.loc[:,'LorBin'].tolist()
        comebin = df.loc[:,'COMEBin'].tolist()
        semibin = df.loc[:,'SemiBin'].tolist()
        semibin2 = df.loc[:,'SemiBin2'].tolist()
        vamb = df.loc[:,'VAMB'].tolist()
        aamb = df.loc[:,'AAMB'].tolist()
        metabat2 = df.loc[:,'MetaBAT2'].tolist()
        y.append(np.mean(lorbin))
        y.append(np.mean(comebin))
        y.append(np.mean(semibin))
        y.append(np.mean(semibin2))
        y.append(np.mean(vamb))
        y.append(np.mean(aamb))
        y.append(np.mean(metabat2))

        y_e.append(np.std(lorbin))
        y_e.append(np.std(comebin))
        y_e.append(np.std(semibin))
        y_e.append(np.std(semibin2))
        y_e.append(np.std(vamb))
        y_e.append(np.std(aamb))
        y_e.append(np.std(metabat2))
        print(x)
        print(y)
        #print(y_e)

        plt.barh(x, y, height=0.6, color=discrete_colors,xerr=y_e,capsize=2,error_kw={'linewidth': 0.5})
        ax.axvline(x=0, color='black', linewidth=1)
        for i, col in enumerate(df.columns):
            # 添加x轴抖动 (jitter)
            y_jittered = np.random.normal(i, 0.02, size=len(df))
            ax.scatter(df[col], y_jittered, color='black', alpha=0.6, s=0.5, zorder=10)

        if n == 0:
            plt.yticks(x, fontsize=8)
        else:
            plt.yticks([])
        plt.xticks([-0.25, 0.25], fontsize=8)
        """
        if env == "airways":
            plt.title("呼吸道", fontsize=12)
        elif env == "oral":
            plt.title("口腔", fontsize=12)
        elif env == "urog":
            plt.title("泌尿生殖道", fontsize=12)
        elif env == "skin":
            plt.title("皮肤", fontsize=12)
        elif env == "gi":
            plt.title("胃肠道", fontsize=12)
            
    """
    plt.figtext(0.38, -0.2,  "Aera under sihouette coefficient curve")
    plt.savefig('areas.svg', bbox_inches='tight', pad_inches=0.02,format='svg')
    plt.show()
