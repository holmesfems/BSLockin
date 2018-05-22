import scipy.interpolate
import scipy.optimize
import numpy
import sys
import os
import re
import datetime

def linear_fmin(func,start,end,step):
    min = func(start)
    xmin = start
    for i in numpy.arange(start,end,step):
        val = func(i)
        if min > val:
            min = val
            xmin = i
        print("i=",i,"func=",val,flush=True) 
    return xmin

def lockin(bsFile, oppsFile, mkidFile,dt = 0.1, mkidDefaultShift = 2.1, BSFreq = 10.0, DAQFreq = 80000.0, FPSFreq = 64.0, BSSE_TimeSpanRate = 10, mkidShiftErr = 0.2,mkidShiftStep = 0.001, ConvolStep = 0.02):
    bsArr = numpy.loadtxt(bsFile)
    oppsArr = numpy.loadtxt(oppsFile)
    mkidArr = numpy.loadtxt(mkidFile)
    mkidArr = numpy.vstack((mkidArr[:,0],mkidArr[:,2])).T
    #Calculate calibrated DAQ Sampling Frequency
    firstPulse = -1
    for i in range(0,len(oppsArr)):
        if (oppsArr[i][1] == 1) and (oppsArr[i+1][1] == 1) and (( -oppsArr[i][0] + oppsArr[i+1][0]) / DAQFreq < 0.1):
            firstPulse = oppsArr[i][0]
            break
    lastPulse = -1
    for i in range(0,len(oppsArr)):
        i = len(oppsArr) - i - 2
        if (oppsArr[i][1] == 1) and (oppsArr[i+1][1] == 1) and (( -oppsArr[i][0] + oppsArr[i+1][0]) / DAQFreq < 0.1):
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
    for i in range(0,len(bsArr_new)):
        if(bsArr_new[i][1] == bsArr_new[i+1][1]) and (-bsArr_new[i][0] + bsArr_new[i+1][0]) * BSFreq < BSSE_TimeSpanRate:
            bs_start = bsArr_new[i][0]
            break
    bs_end = -1.0
    for i in range(0,len(bsArr_new)):
        i = len(bsArr_new) - i - 2
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
    interp_mkid = scipy.interpolate.interp1d(mkidArr[:,0],mkidArr[:,1])
    interp_BS = scipy.interpolate.interp1d(bsArr_new[:,0],bsArr_new[:,1])
    def S(A):
        sumS = 0.0
        for x in numpy.arange(bs_start,bs_end,ConvolStep):
            sumS = sumS + interp_mkid(x - A)*interp_BS(x) * ConvolStep
        return sumS
    mkidShift = linear_fmin(S,mkidDefaultShift,mkidDefaultShift+mkidShiftErr,mkidShiftStep)
    mkidShift_calib = scipy.optimize.fmin(S,mkidShift)
    print("Mkid shift time is:",mkidShift_calib)
    mkidArr = mkidArr + [mkidShift_calib[0],0]
    with open(mkidFile + ".calib","wb") as ofs:
        numpy.savetxt(ofs,mkidArr)
    print("Output calibrated mkid data done!")
    return 0

if not len(sys.argv) < 4:
    lockin(sys.argv[1],sys.argv[2],sys.argv[3])
