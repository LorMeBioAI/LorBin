import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import ScalarFormatter
import matplotlib.scale as mscale
import matplotlib.transforms as mtransforms
import matplotlib.ticker as ticker



class SquareRootScale(mscale.ScaleBase):
    """
    ScaleBase class for generating square root scale.
    """

    name = 'squareroot'

    def __init__(self, axis, **kwargs):
        # note in older versions of matplotlib (<3.1), this worked fine.
        # mscale.ScaleBase.__init__(self)

        # In newer versions (>=3.1), you also need to pass in `axis` as an arg
        mscale.ScaleBase.__init__(self, axis)

    def set_default_locators_and_formatters(self, axis):
        axis.set_major_locator(ticker.AutoLocator())
        axis.set_major_formatter(ticker.ScalarFormatter())
        axis.set_minor_locator(ticker.NullLocator())
        axis.set_minor_formatter(ticker.NullFormatter())

    def limit_range_for_scale(self, vmin, vmax, minpos):
        return max(0., vmin), vmax

    class SquareRootTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def transform_non_affine(self, a):
            return np.array(a) ** 0.5

        def inverted(self):
            return SquareRootScale.InvertedSquareRootTransform()

    class InvertedSquareRootTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def transform(self, a):
            return np.array(a) ** 2

        def inverted(self):
            return SquareRootScale.SquareRootTransform()

    def get_transform(self):
        return self.SquareRootTransform()


mscale.register_scale(SquareRootScale)


if __name__=="__main__":
    basedir = "./realdata/binSize"
    lorbin = (pd.read_csv(f'{basedir}/lorbin.csv', index_col=0).values).reshape(-1)
    semibin = (pd.read_csv(f'{basedir}/semibin.csv', index_col=0).values).reshape(-1)
    semibin2 = (pd.read_csv(f'{basedir}/semibin2.csv', index_col=0).values).reshape(-1)
    aamb = (pd.read_csv(f'{basedir}/aamb_2.csv', index_col=0).values).reshape(-1)
    vamb = (pd.read_csv(f'{basedir}/vamb_2.csv', index_col=0).values).reshape(-1)
    metabat2 = (pd.read_csv(f'{basedir}/metabat2.csv', index_col=0).values).reshape(-1)
    comebin = (pd.read_csv(f'{basedir}/comebin.csv', index_col=0).values).reshape(-1)

    lorbin_total = lorbin.sum()
    semibin_total = semibin.sum()
    semibin2_total = semibin2.sum()
    aamb_total = aamb.sum()
    vamb_total = vamb.sum()
    metabat2_total = metabat2.sum()
    comebin_total = comebin.sum()
    print((lorbin_total-semibin2_total)/semibin2_total*100)
    plt.figure(figsize=(2.4, 2.4), dpi=300)
    x=['LorBin','SemiBin','SemiBin2','VAMB','AAMB','MetaBAT2','COMEBin']
    values =[lorbin_total,semibin_total,semibin2_total,vamb_total,aamb_total,metabat2_total,comebin_total]
    print(values)
    import seaborn as sns

    discrete_colors = sns.color_palette('summer', n_colors=2)
    base_color = discrete_colors[0]
    plt.barh(x, values,color=base_color)
    plt.xscale('squareroot')
    plt.xticks([10,3000,15000,30000],fontsize=8)
    plt.xlabel('# contigs of hmBins (Square root)')
    plt.savefig('contigshmag.png', bbox_inches='tight', pad_inches=0.02)
    plt.show()
"""
# 生成示例数据
data = np.random.randn(1000)

# 绘制平滑的频率分布图
sns.kdeplot(data, shade=True, color='blue')

# 添加标题和标签
plt.title('Smoothed Frequency Distribution')
plt.xlabel('Value')
plt.ylabel('Density')

# 显示图形
plt.show()
"""
