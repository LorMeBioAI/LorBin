import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer
from sklearn.model_selection import train_test_split
from torch.optim.lr_scheduler import StepLR
from torch.utils.data import DataLoader, TensorDataset


class ChannelAttention(nn.Module):
    def __init__(self, in_channels, reduction_ratio=2):
        super(ChannelAttention, self).__init__()
        self.avg_pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Sequential(
            nn.Linear(in_channels, in_channels // reduction_ratio, bias=False),
            nn.ReLU(inplace=True),
            nn.Linear(in_channels // reduction_ratio, in_channels, bias=False), 
        )
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        # 在第二个维度上进行平均池化
        y = self.avg_pool(x.unsqueeze(-1)).squeeze(-1)
        # 应用 fc 层
        y = self.fc(y)
        # 将输出扩展为与输入相同的维度
        return x * self.sigmoid(y)  # 应用 sigmoid 并乘以原始输入 x


class CBAM(nn.Module):
    def __init__(self, in_channels, reduction_ratio=2):
        super(CBAM, self).__init__()
        self.channel_attention = ChannelAttention(in_channels, reduction_ratio)

    def forward(self, x):
        x = self.channel_attention(x)
        return x



class CrossNet(nn.Module):


    def __init__(self, in_features, layer_num=2, seed=1024, device='cpu'):
        super(CrossNet, self).__init__()
        self.layer_num = layer_num
        self.kernels = torch.nn.ParameterList(
            [nn.Parameter(nn.init.xavier_normal_(torch.empty(in_features, 1))) for i in range(self.layer_num)])
        self.bias = torch.nn.ParameterList(
            [nn.Parameter(nn.init.zeros_(torch.empty(in_features, 1))) for i in range(self.layer_num)])

        self.device = device
        self.to(self.device)

    def forward(self, inputs):
        x_0 = inputs.unsqueeze(2)
        x_l = x_0
        for i in range(self.layer_num):
            xl_w = torch.tensordot(x_l, self.kernels[i], dims=([1], [0]))
            dot_ = torch.matmul(x_0, xl_w)
            x_l = dot_ + self.bias[i] + x_l
        x_l = torch.squeeze(x_l, dim=2)
        return x_l

class EvaluationModel(nn.Module):
    def __init__(self, input_size,device='cpu'):
        super(EvaluationModel, self).__init__()
        self.layer=nn.Sequential(
            nn.Linear(input_size*2,input_size*2),
            nn.LeakyReLU(),
            nn.Linear(input_size*2,1),
            nn.Sigmoid(),
        )
        self.device = device
        self.CBAM = CBAM(input_size*2)
        self.CrossNet = CrossNet(input_size,device=self.device)
        self.to(self.device)


    def forward(self, x):
        cross_input = x
        cross_output = self.CrossNet(cross_input)
        stack_out = torch.cat((x, cross_output), dim=-1)
        x = self.CBAM(stack_out)
        return self.layer(x)

            
