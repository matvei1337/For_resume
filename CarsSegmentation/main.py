import os

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
import torch
import torch.utils.data as data
import torchvision
from torchvision import models
import torchvision.transforms.v2 as tfs
import torch.nn as nn
import torch.optim as optim

class SegmentDataset(data.Dataset):
    def __init__(self, path, transform_img=None, transform_mask=None):
        self.path = path
        self.transform_img = transform_img
        self.transform_mask = transform_mask

        path = os.path.join(self.path, "images")
        list_files = os.listdir(path)
        self.length = len(list_files)
        self.images = list(map(lambda _x: os.path.join(path, _x), list_files))
        
        path = os.path.join(self.path, "masks")
        list_files = os.listdir(path)
        self.masks = list(map(lambda _x: os.path.join(path, _x), list_files))
    
    def __getitem__(self, item):
        path_img, path_mask = self.images[item], self.masks[item]
        img = Image.open(path_img).convert("RGB")
        mask = Image.open(path_mask).convert("L")

        if self.transform_img:
            img = self.transform_img(img)
        
        if self.transform_mask:
            mask = self.transform_mask(mask)
            mask[mask < 250] = 1
            mask[mask >= 250] = 0
        
        return img, mask
    
    def __len__(self):
        return self.length

class UnetModel(nn.Module):
    class _TwoConvLayers(nn.Module):
        def __init__(self, in_channels, out_channels):
            super().__init__()
            self.model = nn.Sequential(
                nn.Conv2d(in_channels, out_channels, 3, padding=1, bias=False),
                nn.ReLU(inplace=True),
                nn.BatchNorm2d(out_channels),
                nn.Conv2d(out_channels, out_channels, 3, padding=1, bias=False),
                nn.ReLU(inplace=True),
                nn.BatchNorm2d(out_channels)
            )
        
        def forward(self, x):
            x = self.model(x)
            return x
    
    class _EncoderBlock(nn.Module):
        def __init__(self, in_channels, out_channels):
            super().__init__()
            self.block = UnetModel._TwoConvLayers(in_channels, out_channels)
            self.max_pool = nn.MaxPool2d(2)
        
        def forward(self, x):
            x = self.block(x)
            y = self.max_pool(x)
            return y, x
    
    class _DecoderBlock(nn.Module):
        def __init__(self, in_channels, out_channels):
            super().__init__()
            self.transpose = nn.ConvTranspose2d(in_channels, out_channels, 2, stride=2)
            self.block = UnetModel._TwoConvLayers(in_channels, out_channels)
        
        def forward(self, x, y):
            x = self.transpose(x)
            u = torch.cat([x, y], dim=1)
            u = self.block(u)
            return u



    def __init__(self, in_channels=3, num_classes=1):
        super().__init__()
        self.enc_block1 = self._EncoderBlock(in_channels, 64)
        self.enc_block2 = self._EncoderBlock(64, 128)
        self.enc_block3 = self._EncoderBlock(128, 256)
        self.enc_block4 = self._EncoderBlock(256, 512)

        self.bottleneck = self._TwoConvLayers(512, 1024)

        self.dec_block1 = self._DecoderBlock(1024, 512)
        self.dec_block2 = self._DecoderBlock(512, 256)
        self.dec_block3 = self._DecoderBlock(256, 128)
        self.dec_block4 = self._DecoderBlock(128, 64)

        self.out = nn.Conv2d(64, num_classes, 1)
    
    def forward(self,  x):
        x, y1 = self.enc_block1(x)
        x, y2 = self.enc_block2(x)
        x, y3 = self.enc_block3(x)
        x, y4 = self.enc_block4(x)

        x = self.bottleneck(x)

        x = self.dec_block1(x, y4)
        x = self.dec_block2(x, y3)
        x = self.dec_block3(x, y2)
        x = self.dec_block4(x, y1)

        return self.out(x)


class SoftDiceLoss(nn.Module):
    def __init__(self, smooth=1):
        super().__init__()
        self.smooth = smooth
    
    def forward(self, logits, targets):
        num = targets.shape[0]
        probs = nn.functional.sigmoid(logits)
        m1 = probs.view(num, -1)
        m2 = targets.view(num, -1)
        intersection = (m1 * m2)
        score = 2 * (intersection.sum(1) + self.smooth)/(m1.sum(1) + m2.sum(1) + self.smooth)
        score = 1 - score.sum()/num
        return score


tr_img = tfs.Compose([tfs.ToImage(), tfs.ToDtype(torch.float32, scale=True)])
tr_mask = tfs.Compose([tfs.ToImage(), tfs.ToDtype(torch.float32)])

d_train = SegmentDataset("dataset_seg", transform_img=tr_img, transform_mask=tr_mask)
train_data = data.DataLoader(d_train, batch_size=2, shuffle=True)

model = UnetModel()

optimizer = optim.RMSprop(params=model.parameters(), lr=0.001)
loss1 = nn.BCEWithLogitsLoss()
loss2 = SoftDiceLoss()

epochs = 10
model.train()
epoch_loss = []
for e in range(epochs):
    loss_mean = 0
    lm_count = 0

    for x_train, y_train in tqdm(train_data):
        predict = model(x_train)
        loss = loss1(predict, y_train) + loss2(predict, y_train)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        lm_count += 1
        loss_mean += loss.item()
    
    epoch_loss.append(loss_mean/lm_count)

st = model.state_dict()
torch.save(st, 'model_unet_seg.tar')

img = Image.open("car_1.jpg").convert('RGB')
img = tr_img(img).unsqueeze(0)

p = model(img).squeeze(0)
x = nn.functional.sigmoid(p.permute(1, 2, 0))
x = x.detach().numpy() * 255
x = np.clip(x, 0, 255).astype('uint8')
plt.imshow(x, cmap='gray')
plt.show()