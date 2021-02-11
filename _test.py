import xdf_reader as xr
import PulseDataAdjuster as pda

exp = xr.xdfReader(file_path='exp12(head)/untitled2.xdf',show=False)
exp.findStreamName()
exp.samplingRateVerify()


eeg = exp.getSeriesAndStamps(streamName = 'Apstest',channel=1)
eeg_target = pda.DataAdjuster(timeSeries=eeg['timeSeries'],timeStamps=eeg['timeStamps'],srate=1000,minutes=10)
eeg_target.findpeak(order=400,show=False,diff_max=1.4,diff_min=0.5)


aux = exp.getSeriesAndStamps(streamName = 'Apstest',channel=2)
aux_target = pda.DataAdjuster(timeSeries=aux['timeSeries'],timeStamps=aux['timeStamps'],srate=1000,minutes=10)
aux_target.findpeak(order=400,show=False,diff_max=1.4,diff_min=0.5)

#pda.compareSignal(peakValidation=False,stamp1=eeg_target.getTimeStamp(),stamp2=aux_target.getTimeStamp(),data1=eeg_target.getOutput(),data2=aux_target.getOutput(),peak1=eeg_target.getPeakIndex(),peak2=aux_target.getPeakIndex())

sr1 = exp.getSeriesAndStamps(streamName = 'BioSigPressure',channel=1)
sr1_target = pda.DataAdjuster(timeSeries=sr1['timeSeries'],timeStamps=sr1['timeStamps'],srate=778,minutes=5)
sr1_target.straightOffsetDelete()
#sr1_target.bessel()
sr1_target.FIR(numtaps=701,show=False)
sr1_target.subFilAndOutput()
#sr1_target.FourierTransform()
sr1_target.movemean()
#sr1_target.movemean(show=False)
sr1_target.findpeak_easy(order=300,look=True,show=False,diff_max=1.5,diff_min=0.4)
#pda.compareSignal(peakValidation=True,stamp1=eeg_target.getTimeStamp(),stamp2=sr1_target.getTimeStamp(),data1=eeg_target.getOutput(),data2=sr1_target.getOutput(),peak1=eeg_target.getPeakIndex(),peak2=sr1_target.getPeakIndex())


sr2 = exp.getSeriesAndStamps(streamName = 'BioSigPressure',channel=2)
sr2_target = pda.DataAdjuster(timeSeries=sr2['timeSeries'],timeStamps=sr2['timeStamps'],srate=778,minutes=5)
sr2_target.straightOffsetDelete()
#sr2_target.reverseSignal()
sr2_target.FIR(numtaps=701,show=False)
sr2_target.subFilAndOutput()
sr2_target.movemean()
#sr2_target.movemean()
#sr2_target.FourierTransform()
sr2_target.findpeak_easy(order=300,look=False,show=False)

pda.compareSignal(show=False,data2_shift =0,peakValidation=True,stamp2=sr1_target.getTimeStamp(),stamp1=eeg_target.getTimeStamp(),data2=sr1_target.getOutput(),data1=eeg_target.getOutput(),peak2=sr1_target.getPeakIndex(),peak1=eeg_target.getPeakIndex())
pda.compareSignal(show=False,data2_shift =0,peakValidation=True,stamp2=sr2_target.getTimeStamp(),stamp1=eeg_target.getTimeStamp(),data2=sr2_target.getOutput(),data1=eeg_target.getOutput(),peak2=sr2_target.getPeakIndex(),peak1=eeg_target.getPeakIndex())
#pda.compareSignal(peakValidation=True,stamp1=sr1_target.getTimeStamp(),stamp2=sr2_target.getTimeStamp(),data1=sr1_target.getOutput(),data2=sr2_target.getOutput(),peak1=sr1_target.getPeakIndex(),peak2=sr2_target.getPeakIndex())
# --------------------------------------------------------

ptt = pda.PTTcalculater(data1=sr1_target.getPeakTime(),data2=sr2_target.getPeakTime())
ptt.calclatePTT()
ptt.movemean(show=False)
ptt1=ptt.getPTT()


ptt = pda.PTTcalculater(data1=eeg_target.getPeakTime(),data2=aux_target.getPeakTime())
ptt.calclatePTT()
ptt.movemean(show=False)
ptt2=ptt.getPTT()

cpptt = pda.comparePTT(ptt_target=ptt1,ptt_ref=ptt2,reverse=False,reverse2=True)
cpptt.showGraph()
