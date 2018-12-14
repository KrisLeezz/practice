#-*-coding=utf-8 -*-
#实现markov
#import numpy as np
#def markov():
#	init_array= np.array([0.6,0.2,0.2])
#	transfer_matrix= np.array([[0.2,0.3,0.5],[0.1,0.6,0.3],[0.4,0.5,0.1]])
#	restmp=init_array
#	for i in range(60):
#		res=np.dot(restmp,transfer_matrix)
#		print i, "\t",res
#		restmp=res
#markov()
#傅里叶变换
import numpy as np
from scipy.fftpack import fft,ifft
import matplotlib.pyplot as plt
import seaborn
import csv

def show(ori_func, ft, sampling_period = 5): 
    n = len(ori_func) 
    interval = sampling_period / n 
    # 绘制原始函数
    plt.subplot(2, 1, 1) 
    plt.plot(np.arange(0, sampling_period, interval), ori_func, 'black') 
    plt.xlabel('Time'), plt.ylabel('Amplitude') 
    # 绘制变换后的函数
    plt.subplot(2,1,2) 
    frequency = np.arange(n / 2) / (n * interval) 
    nfft = abs(ft[range(int(n / 2))] / n ) 
    plt.plot(frequency, nfft, 'red') 
    plt.xlabel('Freq (Hz)'), plt.ylabel('Amp. Spectrum') 
    plt.show()


#读取观测数据
X=[]
csv_file_read=open('C:/Users/Mason/Documents/vs2017/new_2015_2017_1.csv')
csv_read=csv.reader(csv_file_read)
for row in csv_read:
	if row[7]!='PM2.5':
		X.append(row[7])
		
x=np.array(X)
x=x.astype(float)

#Y采样数据
y=x

yy=np.fft.fft(y)                     #快速傅里叶变换
yreal = yy.real               # 获取实数部分
yimag = yy.imag               # 获取虚数部分

show(y,yy)
plt.show()

