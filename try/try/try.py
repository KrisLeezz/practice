#-*-coding=utf-8 -*-
#实现markov
import numpy as np
from scipy.fftpack import fft,ifft
import pywt
from PyEMD import EMD, Visualisation
import matplotlib.pyplot as plt
import seaborn
import csv
#马尔可夫方法
def markov():
	init_array= np.array([0.6,0.2,0.2])
	transfer_matrix= np.array([[0.2,0.3,0.5],[0.1,0.6,0.3],[0.4,0.5,0.1]])
	restmp=init_array
	for i in range(60):
		res=np.dot(restmp,transfer_matrix)
		print i, "\t",res
		restmp=res
#markov()

#傅里叶变换
def show(ori_func, ft, sampling_period = 448): 
    n = len(ori_func) 
    interval = sampling_period / n 
    # 绘制原始函数
    plt.subplot(2, 1, 1) 
    plt.plot(np.arange(0, sampling_period, interval), ori_func, 'black') 
    plt.xlabel('Time'), plt.ylabel('RH') 
    # 绘制变换后的函数
    plt.subplot(2,1,2) 
    frequency = np.arange(n / 2) / (n * interval) 
    nfft = abs(ft[range(int(n / 2))] / n ) 
    plt.plot(frequency, nfft, 'maroon')
    #plt.scatter(frequency[0],nfft[0],marker='o',c='maroon',edgecolors='red',linewidths=2)
    plt.xlabel('Freq'), plt.ylabel('RH') 
    plt.show()

#读取观测数据
X=[]
csv_file_read=open('C:/Users/Mason/Documents/vs2017/new_2015_2017_1.csv')
csv_read=csv.reader(csv_file_read)
for row in csv_read:
	if row[13]!='rh':
		X.append(row[13])
		
x=np.array(X)
x=x.astype(float)

#time = np.arange(0, 5, .005) 
#x = np.sin(2 * np.pi * 1 * time)
#x2 = np.sin(2 * np.pi * 2 * time) 
#x3 = np.sin(2 * np.pi * 60 * time) 
#x += x2 +x3
#Y采样数据
y=x
#时间天数
t= np.arange(0,448,1)
def fft(a):
    yy=np.fft.fft(y)#快速傅里叶变换
    return yy
#show(y,yy)
#plt.show()

#小波变换
print pywt.families()
print pywt.wavelist(kind='continuous')
print pywt.ContinuousWavelet('morl')
mode = pywt.Modes.smooth
#信号分解
def plot_signal_decomp(data, w, title):
    """Decompose and plot a signal S.
    S = An + Dn + Dn-1 + ... + D1
    """
    w = pywt.Wavelet(w)#选取小波函数
    a = data
    ca = []#近似分量/低频
    cd = []#细节分量/高频
    for i in range(5):
        (a, d) = pywt.dwt(a, w, mode)#进行5阶离散小波变换
        ca.append(a)
        cd.append(d)
	#展示分量
    rec_a = []
    rec_d = []

    for i, coeff in enumerate(ca):
        coeff_list = [coeff, None] + [None] * i
        rec_a.append(pywt.waverec(coeff_list, w))#反变换

    for i, coeff in enumerate(cd):
        coeff_list = [None, coeff] + [None] * i
        if i ==3:
            print(len(coeff))
            print(len(coeff_list))
        rec_d.append(pywt.waverec(coeff_list, w))
	#作图
    fig = plt.figure()
    ax_main = fig.add_subplot(len(rec_a) + 1, 1, 1)
    ax_main.set_title(title)
    ax_main.plot(data)
    ax_main.set_xlim(0, len(data) - 1)

    for i, y in enumerate(rec_a):
        ax = fig.add_subplot(len(rec_a) + 1, 2, 3 + i * 2)
        ax.plot(y, 'r')
        ax.set_xlim(0, len(y) - 1)
        ax.set_ylabel("A%d" % (i + 1))

    for i, y in enumerate(rec_d):
        ax = fig.add_subplot(len(rec_d) + 1, 2, 4 + i * 2)
        ax.plot(y, 'g')
        ax.set_xlim(0, len(y) - 1)
        ax.set_ylabel("D%d" % (i + 1))
#plot_signal_decomp(y,'rbio2.6','rh - Symmlets5')
#经验模式分解
def EMD_1(a,t):
    emd=EMD()
    emd.emd(a)
    imfs,res=emd.get_imfs_and_residue()
    vis=Visualisation(emd)
    vis.plot_imfs(t=t)
    vis.plot_instant_freq(t)
    vis.show()

IMFS=EMD_1(y,t)
#plt.show()
