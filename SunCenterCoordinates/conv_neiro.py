import os
import json
from PIL import Image

import torch
import torch.utils.data as data
import torchvision.transforms.v2 as tfs
import torch.nn as nn
import torch.optim as optim
from tqdm import tqdm

class SunDataset(data.Dataset):
    def __init__(self, path, train, transform=None):
        self.path = os.path.join(path, "train" if train=="train" else "test")
        self.transform = transform

        with open(os.path.join(self.path, "format.json"), "r") as f:
            self.format = json.load(f)
        
        self.length = len(self.format)
        self.files = tuple(self.format.keys())
        self.targets = tuple(self.format.values())
    
    def __getitem__(self, item):
        path_file = os.path.join(self.path, self.files[item])
        img = Image.open(path_file).convert("RGB")

        if self.transform:
            img = self.transform(img)
        
        return img, torch.tensor(self.targets[item], dtype=torch.float32)
    
    def __len__(self):
        return self.length

transforms = tfs.Compose([tfs.ToImage(), tfs.ToDtype(torch.float32, scale=True)])
d_train = SunDataset("dataset_reg", "train", transforms)
train_data = data.DataLoader(d_train, batch_size=32, shuffle=True)

class SunModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = nn.Conv2d(3, 32, 3, padding='same')#(batch, 3, 256, 256)
        self.pool1 = nn.MaxPool2d(2)
        self.conv2 = nn.Conv2d(32, 8, 3, padding='same')#(batch, 32, 128, 128)
        self.pool2 = nn.MaxPool2d(2)
        self.conv3 = nn.Conv2d(8, 4, 3, padding='same') #(batch, 8, 64, 64)
        self.pool3 = nn.MaxPool2d(2)                    #(batch, 4, 32, 32)
        self.flatten = nn.Flatten()                     #(batch, 4096)
        self.fc1 = nn.Linear(4096, 128)
        self.fc2 = nn.Linear(128, 2)
    
    def forward(self, x):
        x = self.pool1(nn.functional.relu((self.conv1(x))))
        x = self.pool2(nn.functional.relu((self.conv2(x))))
        x = self.pool3(nn.functional.relu((self.conv3(x))))
        x = self.flatten(x)
        x = self.fc2(nn.functional.relu(self.fc1(x)))
        return x

model = SunModel()
optimizer = optim.Adam(params=model.parameters(), lr=0.01, weight_decay=0.001)
loss_func = nn.MSELoss()

epochs = 5
model.train()
losses = []
for e in range(epochs):
    loss_epoch = 0
    lm_count = 0

    for x_train, y_train in tqdm(train_data, leave=True):
        predict = model(x_train)
        loss = loss_func(predict, y_train)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        lm_count += 1
        loss_epoch += loss.item()

    losses.append(loss_epoch/lm_count)

st = model.state_dict()
torch.save(st, "sun_model.tar")

print(losses)

d_test = SunDataset("dataset_reg", "test", transforms)
test_data = data.DataLoader(d_test, batch_size=50, shuffle=True)

Q = 0
count = 0
model.eval()

for x_test, y_test in tqdm(test_data, leave=False):
    with torch.no_grad():
        p = model(x_test)
        Q += loss_func(p, y_test).item()
        count += 1

Q /= count
print(Q)
