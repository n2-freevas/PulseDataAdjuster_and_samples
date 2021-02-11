import xdf_reader as xr
import PulseDataAdjuster as pda

exp = xr.xdfReader(file_path='exp12(head)/untitled4.xdf',show=False)
exp.findStreamName()
exp.samplingRateVerify()

poly = exp.getSeriesAndStamps(streamName = 'Apstest',channel=1)
poly_target = pda.DataAdjuster(timeSeries=poly['timeSeries'],timeStamps=poly['timeStamps'],srate=1000,minutes=3)
poly_target.findpeak(order=500,show=False,diff_max=1.4,diff_min=0.5)


# SRセンサ
sr = exp.getSeriesAndStamps(streamName = 'BioSigPressure',channel=1)
sr_target = pda.DataAdjuster(timeSeries=sr['timeSeries'],timeStamps=sr['timeStamps'],srate=778,minutes=3)
sr_target.straightOffsetDelete()
sr_target.movemean()
#sr_target.bessel()
sr_target.FIR(numtaps=701,show=False)
sr_target.subFilAndOutput()
sr_target.movemean()
#sr_target.showGraph(compare=False)
#sr_target.FourierTransform()
#sr_target.findmaxima(show=False)
sr_target.findpeak_easy(order=300,look=True,show=False,diff_max=1.4,diff_min=0.4)
#sr_target.showGraph(compare=False,appearPeak=True)

# --------------------------------------------------------
pda.compareSignal(data2_shift =0.3,peakValidation=True,stamp2=sr_target.getTimeStamp(),stamp1=poly_target.getTimeStamp(),data2=sr_target.getOutput(),data1=poly_target.getOutput(),peak2=sr_target.getPeakIndex(),peak1=poly_target.getPeakIndex())

