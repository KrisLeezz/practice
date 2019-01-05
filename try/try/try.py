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
all_x=[]
x=[]
csv_file_read=open('C:/Users/Mason/Documents/vs2017/new_2015_2017_1.csv')
csv_read=csv.reader(csv_file_read)
for row in csv_read:
	if row[7]!='PM2.5':
		x.append(row[7])	
x=np.array(x)
x=x.astype(float)

#Y采样数据
y=x
#y= [0,1,22,3,4,5,45,3,77,9,10,11,56,13,14,7,16,17]
#时间天数
t= np.arange(0,len(y),1)
def fft(a):
    yy=np.fft.fft(y)#快速傅里叶变换
    return yy
#show(y,yy)
#plt.show()

#小波变换
print pywt.families()
print pywt.wavelist(kind='discrete')#选择离散还是连续
print pywt.ContinuousWavelet('morl')#连续小波morl的系列
mode = pywt.Modes.smooth
#离散小波变换

def plot_signal_decomp(data, w, title,order):
    """Decompose and plot a signal S.
    S = An + Dn + Dn-1 + ... + D1
    """
    w = pywt.Wavelet(w)#选取小波函数
    a = data
    ca = []#近似分量/低频
    cd = []#细节分量/高频
    for i in range(order):
        (a, d) = pywt.dwt(a, w, mode)#进行order阶离散小波变换
        ca.append(a)
        cd.append(d)
    coeffs = pywt.wavedec(data,w,mode=mode,level=order)# [cA_n, cD_n, cD_n-1, ..., cD2, cD1] : list#与上述方法结果一致
	#各阶重构
    rec_a = []
    rec_d = []

    for i, coeff in enumerate(ca):
        coeff_list = [coeff, None] + [None] * i
        rec_a.append(pywt.waverec(coeff_list, w))

    for i, coeff in enumerate(cd):
        coeff_list = [None, coeff] + [None] * i
        if i ==3:
            print(len(coeff))
            print(len(coeff_list))
        rec_d.append(pywt.waverec(coeff_list, w))
    #slice rec_d
    rec_wave=pywt.waverec(coeffs,w,mode=mode)#重建原始信号=A5+D5+D4+...+D1

	#作图
    fig = plt.figure()
    ax_main = fig.add_subplot(len(rec_a) + 1, 1, 1)
    ax_main.set_title(title)
    ax_main.plot(data,'k-',label='origin wave')
    #ax_main.plot(rec_wave,'c-.',label='reconstruct wave')
    ax_main.legend(loc =1, prop = {'size':5})
    #ax_main.plot(data-rec_wave,'c--')
    ax_main.set_xlim(0, len(data) - 1)

    for i, y in enumerate(rec_a):
        ax = fig.add_subplot(len(rec_a) + 1, 2, 3 + i * 2)
        ax.plot(y, 'r')
        ax.set_xlim(0, len(y) - 1)
        ax.set_ylabel("A%d" % (i + 1),labelpad=1)

    for i, y in enumerate(rec_d):
        ax = fig.add_subplot(len(rec_d) + 1, 2, 4 + i * 2)
        ax.plot(y, 'g')
        ax.set_xlim(0, len(y) - 1)
        ax.set_ylabel("D%d" % (i + 1),labelpad=1)
    #fig.tight_layout()#调整整体空白
    plt.subplots_adjust(wspace =0.3, hspace =1)#调整子图间距
    plt.show()
plot_signal_decomp(y,'db1','PM2.5-dwt',7)
#连续小波变换
def plot_signal_decomp2(data,w,scale,title):
	w = pywt.ContinuousWavelet(w)#选取小波函数
	a = data
	coefs=[]
	freqs=[]
	coefs, freqs=pywt.cwt(a,scale,w)
	#X,Y=coefs.shape#增加等高线表示
	im=plt.matshow(coefs,interpolation='nearest',cmap='jet')
	plt.colorbar(im)
	#plt.title('RH CWT')
	plt.xlabel('time(day)')
	plt.ylabel('scale')
	plt.title(title,verticalalignment='bottom')
	plt.show()
	return coefs
#plot_signal_decomp2(y,'morl',np.arange(1,400),title='RH')
#aod_cwt=[]
#pm_cwt=[]
#so2_cwt=[]
#no2_cwt=[]
#co_cwt=[]
#o3_cwt=[]
#air_cwt=[]
#rh_cwt=[]
#ws_cwt=[]
#p_cwt=[]
#wd_cwt=[]
#x_cwt=[]
#aod_cwt=plot_signal_decomp2(x_t[0],'morl',np.arange(1,400),'AOD_CWT')
#pm_cwt=plot_signal_decomp2(x_t[1],'morl',np.arange(1,400),'PM_CWT')
#so2_cwt=plot_signal_decomp2(x_t[2],'morl',np.arange(1,400),'SO2_CWT')
#no2_cwt=plot_signal_decomp2(x_t[3],'morl',np.arange(1,400),'NO2_CWT')
#co_cwt=plot_signal_decomp2(x_t[4],'morl',np.arange(1,400),'CO_CWT')
#o3_cwt=plot_signal_decomp2(x_t[5],'morl',np.arange(1,400),'O3_CWT')
#air_cwt=plot_signal_decomp2(x_t[6],'morl',np.arange(1,400),'AT_CWT')
#rh_cwt=plot_signal_decomp2(x_t[7],'morl',np.arange(1,400),'RH_CWT')
#ws_cwt=plot_signal_decomp2(x_t[8],'morl',np.arange(1,400),'WS_CWT')
#p_cwt=plot_signal_decomp2(x_t[9],'morl',np.arange(1,400),'P_CWT')
#wd_cwt=plot_signal_decomp2(x_t[10],'morl',np.arange(1,400),'WD_CWT')
#x_cwt.append(aod_cwt)
#x_cwt.append(pm_cwt)

#经验模式分解
def EMD_1(a,t):
    emd=EMD()
    emd.emd(a)
    imfs,res=emd.get_imfs_and_residue()
    plt.plot(t,a)
    plt.title('origin sequence')
    vis=Visualisation(emd)
    vis.plot_imfs(t=t)
    vis.plot_instant_freq(t)
    vis.show()
    plt.show()
#IMFS=EMD_1(y,t)
#plt.show()
