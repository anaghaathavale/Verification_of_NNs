arr = list(map(float, input().split()))

import numpy as np

#x1 = input()
#x2 = input()
#x3 = input()
#x4 = input()

inputs = arr

n11_wts = [-0.71879154,0.10777157, 0.61765987, 0.32582593]
n12_wts = [0.06999272,1.223484,-0.6634537, -0.65185314] 
n13_wts = [-0.5479262,0.30571675,-0.427665,0.06441003]
b11 = -0.08128614,
b12 = 0.44203088
b13 = 0

n21_wts =[-0.31779695, -0.45176014, -0.3360595]
n22_wts =[0.6303168, 1.0507578 ,0.5606068]
n23_wts = [0.24277228, 1.2744203 ,-0.22409433]
n24_wts = [-0.34698159,-0.17203873,0.33977765]
n25_wts = [0.12579483, 1.8299768,0.03665435]
b21 = 0
b22 = -0.8150968
b23 = -0.00104406
b24 = 0
b25 = -0.00138543

n31_wts = [0.67072374, 0.21325482, 0.8423685, 0.3852225, 0.5454857]
n32_wts = [0.38538057,-1.284013,0.7352335,0.36256438,0.6089774]
n33_wts = [0.00611049,-0.8554969,-0.7638753,-0.41221866,-1.3374699]
b31 = -1.123881
b32 = 0.13953435
b33 = 0.97535825
#print("Input values = ", inputs)
n11 = np.maximum(((np.dot(inputs, n11_wts)) + b11), 0)
n12 = np.maximum(((np.dot(inputs, n12_wts)) + b12), 0)
n13 = np.maximum(((np.dot(inputs, n13_wts)) + b13), 0)
layer1 = [n11, n12, n13]
print("Layer 1 neuron values = ", layer1)

n21 = np.maximum(((np.dot(layer1, n21_wts)) + b21),0)
n22 = np.maximum(((np.dot(layer1, n22_wts)) + b22),0)
n23 = np.maximum(((np.dot(layer1, n23_wts)) + b23),0)
n24 = np.maximum(((np.dot(layer1, n24_wts)) + b24),0)
n25 = np.maximum(((np.dot(layer1, n25_wts)) + b25),0)
layer2 = [n21[0], n22[0], n23[0], n24[0], n25[0]]
print("Layer 2 neuron values = ", layer2)

n31 = np.maximum(((np.dot(layer2, n31_wts)) + b31),0)
n32 = np.maximum(((np.dot(layer2, n32_wts)) + b32),0)
n33 = np.maximum(((np.dot(layer2, n33_wts)) + b33),0)
layer3 = [n31, n32, n33]
print("Layer 3 neuron values = ", layer3)
