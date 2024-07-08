import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from torch.optim.lr_scheduler import StepLR
from torch.utils.data import DataLoader, TensorDataset



class KeepModel(nn.Module):
    def __init__(self, input_size):
        super(KeepModel, self).__init__()
        self.linear1 = nn.Linear(input_size,input_size)
        self.relu = nn.LeakyReLU()
        self.linear2 = nn.Linear(input_size,1)
        self.sigmoid = nn.Sigmoid()
    def forward(self, x):
        x = self.relu(self.linear1(x))
        x = self.sigmoid(self.linear2(x))
        return x


