import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, fftpack
import warnings
import sys
import csv
import os.path
import codecs
import emd
warnings.filterwarnings("ignore", category=RuntimeWarning)


class DataAdjuster(object):

  def __init__(self,timeSeries=None,timeStamps=None,srate=0,minutes=10):
    
    # 引数をクラス内変数に代入する
    self.datalen = len(timeSeries)    #データのインデックス長さ
    tSe = []
    tSt = []
    for i in range(self.datalen):
      tSe.append(timeSeries[i] + (i*0.00001))
    for i in range(self.datalen):
      tSt.append(timeStamps[i]- timeStamps[0])
    self.series = np.array(tSe)
    self.stamps = np.array(tSt)
    
    self.Fsamp = srate
    self.Tsamp = 1/self.Fsamp
    self.Fnyq = self.Fsamp/2    #ナイキスト周波数
    self.second = minutes*60 #データの計測時間．

    self.average_of_series = np.average(self.series) # series の平均値
    self.series_ave = self.series - self.average_of_series #seriesをseriesの平均値で引いたもの
    self.output = self.series_ave
    self.message('------ import data -------')
    print(self.series_ave)
    print(self.stamps)
    print(self.Fsamp,'Hz /',self.datalen/self.Fsamp,'sec data')
    self.message('----------------------------')
    #各フィルタのパラメータ初期値(呼吸(0.2~0.3Hz)取得用の値に設定)
    self.Fpass = 0.7 #通過域端の周波数  # 1.0
    self.Fstop = 0.5 #阻止域端の周波数 # 0.9
    self.gpass = 0.001 #通過域最大損失量[dB]
    self.gstop = 0.05 #阻止域最小減衰量[dB]

    #正規化
    self.Wpass = self.Fpass/self.Fnyq
    self.Wstop = self.Fstop/self.Fnyq
    self.peakindex = []
    self.peaktime = []
    self.short_peakindex = [] #ピーク間隔が狭すぎる2つのピークの右側のインデックス
    self.short_peaktime = []
    self.long_peakindex = []  #ピーク間隔が広すぎる2つのピークの右側のインデックス
    self.long_peaktime = []

    self.delay = 0
    self.argsCheck()

  # data secures function
  def argsCheck(self):
    if self.datalen != len(self.stamps):
      self.message('[!] The data lengths of Series and Stamps are different.')
    
  #inside system function
  def message(self,message):
    print('[Pulse Data Adjuster] : ',message)

  def straightOffsetDelete(self,show=False,save=True):
    y = self.series_ave
    start = y[0]
    end = y[-1]
    tilt = (end - start)/(self.datalen)
    
    for i in range(self.datalen):
       y[i] -= start + (i * tilt)

    if show:
      self.showGraph(data=y)
    if save:
      self.series_ave = y
      self.output = y
      
  
  def reverseSignal(self,show=False,save=True):
    y = self.series_ave
    y = -1*y
    if show:
      self.showGraph(data=y)
    if save:
      self.series_ave = y
      self.output = y

  # return inclass-variant fuction
  def getPeakTime(self):
    return self.peaktime
  def getPeakIndex(self):
    return self.peakindex

  def getTimeStamp(self):
    return self.stamps
  def getOutput(self):
    return self.output

  #interface function
  def showOutput(self):
    print(self.output)

  def showGraph(self,data=[],appearPeak=False,compare=True):
    self.message('@@ show result')      
    
    if len(data)==0:
      result = self.output
    else:
      result = data
    if compare:
      plt.plot(self.stamps,self.series_ave, label='input')
    plt.plot(self.stamps,result,label='output')
    plt.xlabel('Time[s]')
    plt.ylabel('Ampltude')
    if appearPeak and (len(self.peakindex)!=0):
      plt.plot(self.peaktime,self.output[self.peakindex],'ro',markersize=8,label='peaks(high)/output')
      plt.plot(self.short_peaktime,self.output[self.short_peakindex],'bd',label='peaks(high)/output')
      plt.plot(self.long_peaktime,self.output[self.long_peakindex],'gs',label='peaks(high)/output')
      
    plt.show()

  def FourierTransform(self):        
    F = np.fft.fft(self.output, n=None, axis= -1, norm=None)
    freq = np.fft.fftfreq(self.datalen, d = self.Tsamp)
    Amp = np.abs(F/(self.datalen/2)) # 振幅計算
    
    fig, ax = plt.subplots()
    ax.plot(freq[1:int(self.datalen/2)], Amp[1:int(self.datalen/2)])
    ax.set_xlabel("Frequency [Hz]")
    ax.set_xlim([0,2.0])
    ax.set_ylabel("amplitude")
    ax.grid()
    plt.show()

  def subFilAndOutput(self,show=False,save=True):
    y = self.series_ave - self.output
    if show:
      self.showGraph(data = y)
    if save:
      self.output = y

  def butter(self,save=True,show=False):
    N, Wn = signal.buttord(self.Wpass, self.Wstop, self.gpass,self.gstop)
    b, a = signal.butter(N, Wn, 'low')
    y = signal.filtfilt(b,a,self.output)
    
    if show:
      self.showGraph(data=y)
    if save:
      self.output = y

  def cheby_1(self,save=True,show=False):
    N, Wn = signal.cheb1ord(self.Wpass, self.Wstop, self.gpass,self.gstop)
    b, a = signal.cheby1(N, Wn, 'low')
    y = signal.filtfilt(b,a,self.output)
    
    if show:
      self.showGraph(data=y)
    if save:
      self.output = y
  
  def cheby_2(self,save=True,show=False):
    N, Wn = signal.cheb2ord(self.Wpass, self.Wstop, self.gpass,self.gstop)
    b, a = signal.cheby2(N, Wn, 'low')
    y = signal.filtfilt(b,a,self.output)
    
    if show:
      self.showGraph(data=y)
    if save:
      self.output = y

  def bessel(self,save=True,show=False):
    N, Wn = signal.ellipord(self.Wpass, self.Wstop, self.gpass,self.gstop)
    b, a = signal.bessel(N, Wn, 'low')
    y = signal.filtfilt(b,a,self.output)
    
    if show:
      self.showGraph(data=y)
    if save:
      self.output = y

  def FIR(self,numtaps=1001,show=False,save=True,row_pass=True):
    a = 1
    _numtaps = numtaps
    b = signal.firwin(_numtaps, self.Wpass, window='hann', pass_zero=row_pass,fs=self.Fsamp)
    y = signal.filtfilt(b,a,self.output)
    self.message('@ FIR filtered :'+str(len(y)))

    if show:
      self.showGraph(data=y)
    if save:
      self.output = y
      self.delay = (numtaps-1)/2*self.Tsamp


  def movemean(self,show=False,save=True):
    num = 100
    b = np.ones(num)/num
    y = np.convolve(self.output, b, mode = 'same')

    if show:
      self.showGraph(data=y)
    if save:
      self.output = y
      

  def findmaxima(self,show= False,look=False):
    self.peakindex = signal.argrelmax(self.output,order=1)[0]
    for i in self.peakindex:
      self.peaktime.append(self.stamps[i])
    if show:
      self.showGraph(appearPeak=True)


  def findpeak(self,diff_min=0.6,diff_max=1.5,order=325,show = False,look=False):
  
    self.peakindex = signal.argrelmax(self.output,order=order)[0]
    for i in self.peakindex:
      self.peaktime.append(self.stamps[i])
    
    #peak_diff = np.diff(self.peaktime,n=1).tolist()
    
    _pass = False
    array = []
    for i in range(len(self.peakindex) - 1):
      if _pass:
        _pass = False
      else:
        array.append(self.peakindex[i])
        now = self.peaktime[i]
        nex = self.peaktime[i+1]
        diff = nex - now
        
        # long peaks を探す
        # 1.5s => 0.666回/1s　=>39ppm →死人？？

        if  diff> diff_max:
          self.long_peakindex.append(self.peakindex[i])
          self.long_peakindex.append(self.peakindex[i+1])
          self.long_peaktime.append(now)
          self.long_peaktime.append(nex)
          array.append(int((self.peakindex[i] + self.peakindex[i+1]) / 2 ))

        # short peaks を探す
        elif diff < diff_min:
          self.short_peakindex.append(self.peakindex[i])
          self.short_peakindex.append(self.peakindex[i+1])
          self.short_peaktime.append(now)
          self.short_peaktime.append(nex)
          _pass = True
        #正常なpeaksの場合
        
    array.append(self.peakindex[-1])

    if print:
      self.message('[find peak] Finish Peaks and error-estimated peaks')
      self.message('pure peaks amount of ' + str(len(self.peakindex)))
      self.message('short pair peaks amount of ' + str(len(self.short_peakindex)))
      self.message('long pair peaks amount of ' + str(len(self.long_peakindex)))
      self.message('finally peaks amount of ' + str(len(array)))

    self.peakindex = np.array(array)
    self.peaktime = []
    for i in self.peakindex:
      self.peaktime.append(self.stamps[i])
    
    if show:
      self.showGraph(appearPeak=True)

  def findpeak_easy(self,diff_min=0.6,diff_max=1.5,order=325,show = False,look=False):
    
    self.peakindex = signal.argrelmax(self.output,order=order)[0]
    for i in self.peakindex:
      self.peaktime.append(self.stamps[i])
    
    if print:
      self.message('[find peak] Finish Peaks and error-estimated peaks')
      self.message('finally peaks amount of ' + str(len(self.peakindex)))

    
    if show:
      self.showGraph(appearPeak=True)


class PTTcalculater(object):
  def __init__(self,data1=None, data2=None):
    if(data1==None) or (data2==None):
      self.message('Please prepare 2 data array of peaks')
      return None
    
    # data1 が先に発生する(心臓に近い)
    self.ptt_pre = []
    self.ptt = []
    self.pttstamp = []
    self.data1 = data1
    self.data2 = data2
    self.data1len = len(self.data1)
    self.data2len = len(self.data2)
    self.datalen = 0
    self.message('data1 is:'+str(self.data1len))
    self.message('data2 is:'+str(self.data2len))

    #　最長のデータ長をもつ配列のそれを獲得する
    
  def showData1(self):
    for i in range(self.data1len):
      print(i,self.data1[i])
  def showData2(self):
    for i in range(self.data2len):
      print(i,self.data2[i])

  def checkindex(self,index1,index2):
    if (index1 == self.data1len)or(index2 == self.data2len):
      return False
    else:
      return True

  def calclatePTT(self):
    ptt = []
    pttstamp = []
    index2 = 0
    ptt_ave = 0
    for i in range(self.data1len - 1):
      process = 0
      ptt_candi = []
      while True:
        now = index2 + process
        if now >= self.data2len:
          #print('aaa')
          break
        target = self.data2[now]
        
        if (self.data1[i]<=target)and(target<=self.data1[i+1]):
          
          #print(' /[in]',now,i,target,self.data1[i],target - self.data1[i],end='')
          ptt_candi.append(target - self.data1[i])
          process += 1
        elif(target < self.data1[i]):
          #print(' /[<]',now,i,target,self.data1[i],target - self.data1[i],end='')
          #print(' /process')
          process += 1 
        elif(self.data1[i+1] < target):
          #print(' /[>]',now,i,target,end='')
          #print(' /next :',end='')
          index2 += process
          break
        else:
          print('fuck!!!!!!!!!!!!!!!!!')

      candilen = len(ptt_candi)
      if candilen == 0:
        if len(ptt) > 1:
          #print('inherit',ptt[len(ptt)-1])
          #print('--')
          ptt.append(ptt[len(ptt)-1])
          pttstamp.append(-1)
        else:
          pass
      elif candilen > 1:
        if len(ptt) > 1:
          #print('inherit',ptt[len(ptt)-1])
          #print('--')
          ptt.append(ptt[len(ptt)-1])
          pttstamp.append(-1)
        else:
          pass
      else:
        #print(ptt_candi[0])
        ptt.append(ptt_candi[0])
        pttstamp.append(self.data1[i])
        #print('--')
    
    self.ptt = np.array(ptt)
    
    pttstamp_ave_candiarray = []
    pttstamp_len = len(pttstamp)
    for j in range(1,pttstamp_len-2):
      if(pttstamp[j] == -1):
        if (pttstamp[j-1] != -1) and (pttstamp[j+1] != -1):
          pttstamp[j] = (pttstamp[j-1]+pttstamp[j+1])/2
      else:
        if(pttstamp[j-1] != 1):
          pttstamp_ave_candiarray.append(pttstamp[j] - pttstamp[j-1])
    pttstamp_ave = 0
    for i in pttstamp_ave_candiarray:
      pttstamp_ave += i
    pttstamp_ave = pttstamp_ave/len(pttstamp_ave_candiarray)

    for j in range(1,pttstamp_len):
      if pttstamp[j] == -1:
       if pttstamp[j-1] != -1:
        pttstamp[j] = pttstamp[j-1] + pttstamp_ave
    for j in range(pttstamp_len-1,0,-1):
      if pttstamp[j] == -1:
        if pttstamp[j+1] != -1:
          pttstamp[j] = pttstamp[j+1] - pttstamp_ave

    self.pttstamp = np.array(pttstamp)
    self.datalen = len(self.ptt)

  def EMD(self):
    if len(self.ptt) == 0:
      print('Please calcrates EMD after PTT calcrate')
    else:
      print(type(self.ptt),len(self.ptt))
      imf = emd.sift.sift(self.ptt)
      
      plt.figure()
      print(len(imf))
      horizontal_axis = range(0,len(imf[:,0]))
      plt.plot(horizontal_axis,imf[:,0],'r',label='EMD')
      horizontal_axis = range(0,self.datalen)
      plt.plot(horizontal_axis,self.ptt,'b',label='ptt')
      plt.xlabel('index')
      plt.ylabel('Time[s]')
      plt.show()
  
  def showGraph(self,compare=False):
    plt.figure()
    horizontal_axis = range(0,self.datalen)
    plt.plot(self.stamp1,self.ptt,label='ptt')
    plt.xlabel('index')
    plt.ylabel('Time[s]')
    plt.show()

  def getPTT(self):
    dicti = {'ptt':self.ptt,'stamp':self.pttstamp}
    return dicti

  def movemean(self,show=False,save=True):
    num = 91
    b = np.ones(num)/num
    y = np.convolve(self.ptt, b, mode = 'same')
    y = y[int(num/2)+1:int(num/2)*(-1)]
    self.pttstamp = self.pttstamp[int(num/2)+1:int(num/2)*(-1)]
    leny = len(y)
    if show:
      plt.figure()
      horizontal_axis = range(0,leny)
      plt.plot(horizontal_axis,y,label='ptt')
      plt.xlabel('index')
      plt.ylabel('Time[s]')
      plt.show()
    if save:
      self.ptt = y
      self.datalen = leny
  
  def message(self,message):
    print('[calculatePTT] : ',message)


class comparePTT(object):
  def __init__(self,ptt_ref=None,ptt_target=None,reverse=False,reverse2=False):
    self.ptt1 = np.array(ptt_ref['ptt'])
    self.ptt2 = np.array(ptt_target['ptt'])
    self.stamp1 = np.array(ptt_ref['stamp'])
    self.stamp2 = np.array(ptt_target['stamp'])
    self.ptt1len = len(self.ptt1)
    self.ptt2len = len(self.ptt2)  
    print('[comparePTT] input > ',len(self.ptt1),len(self.ptt2),len(self.stamp1),len(self.stamp2))
    
    if self.ptt1len > self.ptt2len:
      self.ptt1 = self.ptt1[0:self.ptt2len]
      self.ptt1len = self.ptt2len
      self.stamp1 = self.stamp1[0:self.ptt2len]

    elif self.ptt1len < self.ptt2len:
      self.ptt2 = self.ptt2[0:self.ptt1len]
      self.ptt2len = self.ptt1len
      self.stamp2 = self.stamp2[0:self.ptt2len]
    """
    ptt_target_candi = []
    for i in self.stamp1:  
      idx = np.abs(np.asarray(self.stamp2) - i).argmin()
      print(str(idx)+' / ',end='')
      ptt_target_candi.append(self.ptt2[idx])
    self.ptt2 = np.array(ptt_target_candi)
    """
    print('[comparePTT] after corection > ',len(self.ptt1),len(self.ptt2))
    self.ptt1ave = self.ptt1
    self.ptt2ave = self.ptt2
    self.ptt1ave -= self.ptt1ave.mean()
    self.ptt2ave -= self.ptt2ave.mean()   
    if reverse:
      self.ptt2 = -1 * self.ptt2
      self.ptt2ave = -1 * self.ptt2ave 

    if reverse2:
      self.ptt2 = np.flipud(self.ptt2)
      self.ptt2ave = np.flipud(self.ptt2ave)

    print('var:',np.var(self.ptt1ave-self.ptt2ave))
    print('corrcoef:',np.corrcoef(self.ptt1,self.ptt2))
  
  def showGraph(self):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    horizontal_axis = range(0,self.ptt1len)
    ln1 = ax1.plot(horizontal_axis,self.ptt1ave,color='mediumvioletred',label='ECG-pulse')
    ax2 = ax1.twinx()
    ln2 = ax2.plot(horizontal_axis,self.ptt2ave,color='darkblue',label='SRsensor')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc='lower right')
    ax1.set_xlabel('index')
    ax1.set_ylabel('ECG-pulse')
    ax1.grid(True)
    ax2.set_ylabel('SRsensor')
    plt.show()

class compareSignal(object):
  def __init__(self,show=True,stamp1=None,stamp2=None,data1=None,data2=None,peak1=None,peak2=None,peakValidation=False,data2_shift=0):
    
    self.stamp1 = np.array(stamp1)
    self.stamp2 = np.array(stamp2)
    self.data1 = np.array(data1)
    self.data2 = np.array(data2)
    self.data1 -= self.data1.mean()
    self.data2 -= self.data2.mean()
    self.peak1 = np.array(peak1)
    self.peak2 = np.array(peak2)

    peakave1 = self.data1[peak1].mean()
    peakave2 = self.data2[peak2].mean()
    print(peakave1,peakave2)
    
    self.stamp2 = self.stamp2 + data2_shift

    if(peakave1>peakave2):
      _max = peakave1/peakave2
      for i in range(len(self.data2)):
        self.data2[i] = self.data2[i]*_max
    else:
      _max = peakave2/peakave1
      for i in range(len(self.data1)):
          self.data1[i] = self.data1[i]*_max
    
    
    flag = False
    if ((peak1.all())and(peak2.all())):
      flag = True
    
    if show:
      plt.figure()
      plt.plot(self.stamp2,self.data2,label='data2',linestyle='solid',color='mediumblue',linewidth=1)
      plt.xlabel('horizon')
      plt.ylabel('vertical')
      plt.plot(self.stamp1,self.data1,label='data1',linestyle='dashed',color='red',linewidth=1)
      
      if flag:
        plt.plot(self.stamp2[peak2],self.data2[peak2],'bo',alpha=0.7,markersize=6,label='peaks(high)/output')
        plt.plot(self.stamp1[peak1],self.data1[peak1],'ro',alpha=0.7,markersize=6,label='peaks(high)/output')
        
    
    if peakValidation:
      print('[peakValidation]====Start=====')
      clear = 0
      lack = 0
      over = 0
      peak1len = len(self.peak1)
      index = 0
      ptt = []
      for i in range(len(self.peak2)-1):
        que = []
        now = self.stamp2[self.peak2[i]]
        nex = self.stamp2[self.peak2[i+1]]
        p_first = self.stamp1[self.peak1[index]]
        for j in range(index,peak1len):
          p = self.stamp1[self.peak1[j]]
          if (now < p)and(p < nex):
            
            que.append(p)
          elif(nex < p):
            index = j
            break
        if len(que) == 1:
          clear += 1
          ptt.append(p_first - now)
        elif len(que) > 1:
          print('over:',now)
          over += 1
        else:
          print('lack:',now)
          lack += 1

      print('[peakValidation]==========')
      print('clear:',clear)
      print('lack:',lack)
      print('over:',over)
      npppt = np.array(ptt)
      print('ptt_ave:',npppt.mean())
      print('ptt_ave:',np.std(npppt))
      print('==========================')
      
    plt.show()