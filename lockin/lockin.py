import scipy
import numpy
import sys
import os
import re
import datetime

def lockin(bsFile, oppsFile, mkidFile,dt = 0.1, mkidDefaultShift = 0, BSFreq = 10.0, DAQFreq = 80000.0, FPSFreq = 64.0, BSSE_TimeSpanRate = 10):
    bsArr = numpy.loadtxt(bsFile)
    oppsArr = numpy.loadtxt(oppsFile)
    mkidArr = numpy.loadtxt(mkidFile)
    mkidArr = numpy.c_(mkidArr[:,0],mkidArr[:,2])
    #Calculate calibrated DAQ Sampling Frequency
    firstPulse = -1
    for i in range(0,oppsArr.count()):
        if (oppsArr[i][1] == 1) and (oppsArr[i+1][1] == 1) and (( -oppsArr[i][0] + oppsArr[i+1][0]) / DaqFreq < 0.1):
            firstPulse = oppsArr[i][0]
            break
    lastPulse = -1
    for i in range(0,oppsArr.count()):
        i = oppsArr.count() - i - 2
        if (oppsArr[i][1] == 1) and (oppsArr[i+1][1] == 1) and (( -oppsArr[i][0] + oppsArr[i+1][0]) / DaqFreq < 0.1):
            lastPulse = oppsArr[i][0]
            break
    if firstPulse < 0 or lastPulse < 0:
        print("Failed to fit opps pulse!")
        return
    span = round((lastPulse-firstPulse)/DAQFreq)
    DAQFreq_calib = (lastPulse - firstPulse) / span
    
    #Calculate Zero point of time:
    #Zero Point is firstPulse of opps
    #Get first timestamp of BS
    bsFileOpen = open(bsFile);
    timeStr = bsFileOpen.readline()
    timeStrRe = re.compile(r"#(?P<year>\d+)-(?P<month>\d+)-(?P<day>\d+)T(?P<hour>\d+):(?P<minute>\d+):(?P<second>\d+)\.(?P<microsecond>\d+)")
    timeStrMatch = timeStrRe.match(timeStr)
    baseDT = datetime.datetime(int(timeStrMatch.group('year')),\
        int(timeStrMatch.group('month')),\
        int(timeStrMatch.group('day')),\
        int(timeStrMatch.group('hour')),\
        int(timeStrMatch.group('minute')),\
        int(timeStrMatch.group('second')),\
        int(timeStrMatch.group('microsecond')),\
        )
    microTime = float(timeStrMatch.group('microsecond')) / 1.0e6
    zeroPoint = round(microTime + firstPulse*1.0 / DAQFreq_calib)
    timedelta = datetime.timedelta(microseconds = int(round((zeroPoint-microTime)*1.0e6)))
    baseDT_calib = baseDT + timedelta
    #write out new BS array
    bsArr_new = (bsArr - [firstPulse,0])/[DAQFreq_calib,1]
    bsArr_new_file = bsFile + ".calib"
    with open(bsArr_new_file,"w") as ofs:
        ofs.write(baseDT_calib.strftime("#%Y-%m-%dT%H:%M:%S.%f\n"))
    with open(bsArr_new_file,"ab") as ofs:
        numpy.savetxt(ofs,bsArr_new)
    print("Output calibrated beamswitch signal done!")

    #Get time that BS starts
    bs_start = -1.0;
    for i in range(0,bsArr_new.count()):
        if(bsArr_new[i][1] == bsArr_new[i+1][1]) and (-bsArr_new[i][0] + bsArr_new[i+1][0]) * BSFreq < BSSE_TimeSpanRate:
            bs_start = bsArr_new[i][0]
            break
    bs_end = -1.0
    for i in range(0,bsArr_new.count()):
        i = bsArr_new.count() - i - 2
        if(bsArr_new[i][1] == bsArr_new[i+1][1]) and (-bsArr_new[i][0] + bsArr_new[i+1][0]) * BSFreq < BSSE_TimeSpanRate:
            bs_end = bsArr_new[i][0]
            break
    if bs_end > bs_start:
        print("Found BS start point and end point")
        print("BS_Start = ", bs_start)
        print("BS_End = ", bs_end)
    else:
        print("Can't find BS start point and end point")
        return
    #Get calibrated mkid shift
    mkidArr = mkidArr * [1.0/FPSFreq,1.0]
    interp_BS = scipy.interpolate.interp1d(bsArr_new[:,0],bsArr_new[:,1])
    def S(A):
        sumS = 0.0
        sumN = 0
        for item in filter(lambda x: x[0] + A > bs_start and x[0] + A < bs_end, mkidArr):
            sumS = sumS + interp_BS(item[0])*item[1]
            sumN = sumN + 1
        return sumS/sumN
    mkidShift = scipy.optimize.fmin(S,mkidDefaultShift)
    print("Mkid shift time is:",mkidShift)
    mkidArr = mkidArr + [mkidShift,0]
    with open(mkidFile + ".calib","wb") as ofs:
        numpy.savetxt(ofs,mkidArr)
    print("Output calibrated mkid data done!")
    return 0

if not len(sys.argv) < 4:
    lockin(sys.argv[1],sys.argv[2],sys.argv[3])