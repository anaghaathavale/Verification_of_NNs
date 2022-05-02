import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from tensorflow.keras.utils import to_categorical

dataset=pd.read_csv('~/Downloads/Iris.csv')
dataset

import seaborn as sb
sb.set(style="ticks")
sb.set_palette("husl")
sb.pairplot(dataset.iloc[:,1:6],hue="Species")

X = dataset.iloc[:,1:5].values
y = dataset.iloc[:,5].values

from sklearn.preprocessing import LabelEncoder
encoder =  LabelEncoder()
y1 = encoder.fit_transform(y)

Y = pd.get_dummies(y1).values

from sklearn.model_selection import train_test_split  

from keras.models import Sequential
from keras.layers import Dense
model = Sequential()
model.add(Dense(3,input_dim=4,activation='relu'))
model.add(Dense(5,activation='relu'))
model.add(Dense(3,activation='softmax'))

model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])model.summary()

print(model.get_weights())

X_train,X_test, y_train,y_test = train_test_split(X,Y,test_size=0.2,random_state=0) 

model.fit(X_train,y_train,epochs=100,batch_size=8)

model.summary()
print(model.get_weights())

scores = model.evaluate(X_test,y_test)
print(scores[1]*100)
