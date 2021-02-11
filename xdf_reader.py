import pyxdf
import numpy as np


class xdfReader(object):
  def __init__(self,file_path = None,show=False):
    
    self.datas, self.header = pyxdf.load_xdf(file_path)
    if show:
      self.showData()

  def showData(self):
    print('\n[xdfRead]=======showData============\n')
    print('  ∟-> keys name is',self.datas[0].keys())
    for data in self.datas:
      print(data['info'])
      print('series length : ',(len(data['time_series'])))
      print('stamps length : ',(len(data['time_stamps'])))
      print('\n============================\n')

  def samplingRateVerify(self):
    for data in self.datas:
      
      print('[xdfRead][Sampling Rate Verify] : '+data['info']['name'][0]+'\'s data is verified.')
      
      # Information Collect
      stime = round(1/int(data['info']['nominal_srate'][0]),10) # between sampling times
      deviation = 0
      length = len(data['time_stamps'])
      for i in range(length-1):
        between = data['time_stamps'][i+1] - data['time_stamps'][i]
        deviation += abs(between - stime)

      # Varidation Section
      deviation = deviation/length
      print(' ∟-> deviation to nominal sampling rate :'+str(round(deviation,len(str(stime))+3))+' (dT ='+'{:e}'.format(stime)+')\n')
  
  def findStreamName(self,name=None):
    if name==None:
      print('[xdfRead][:)] Streams in this file is there.')
      for i in range(len(self.datas)):
        print(' ∟-> ',self.datas[i]['info']['name'][0],' (Fs=',self.datas[i]['info']['nominal_srate'][0]+') (Chs=',self.datas[i]['info']['channel_count'][0]+')')
    else:
      index = None
      for i in range(len(self.datas)):
        
        if self.datas[i]['info']['name'][0] == name:
          index = i
      if index != None:
        return index
      else:
        ('[xdfRead][!] Stream name you choose not find.')

  def getFs(self,streamName):
    return self.datas[self.findStreamName(streamName)]['info']['nominal_srate']

  def getTimeSeries(self,streamName,channel=None):
    if (channel != None) and (type(channel) is int):
      ch_series = []
      for sample in self.datas[self.findStreamName(streamName)]['time_series']:
        
        ch_series.append(sample[channel - 1])
      returner = np.array(ch_series)
    else:
      returner = self.datas[self.findStreamName(streamName)]['time_series']
    return returner

  def getTimeStamps(self,streamName):
    return self.datas[self.findStreamName(streamName)]['time_stamps']

  def getSeriesAndStamps(self,streamName,channel=None):
    return {'timeSeries':self.getTimeSeries(streamName=streamName,channel=channel),'timeStamps':self.datas[self.findStreamName(streamName)]['time_stamps']}