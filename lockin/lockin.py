import scipy.interpolate
import scipy.optimize
import scipy.integrate
import scipy.signal
import numpy
import sys
import os
import re
import datetime
import math
import progbar
import time
#Params
lockin_dt = 0.25
mkidDefaultShift = 0.0
BSFreq = 10.0
DAQFreq = 80000.0
FPSFreq = 64.0
BSSE_TimeSpanRate = 10
LSE_TimeSpanRate = 0.55
mkidShiftErr = 0.5
mkidShiftStep = 0.001
ConvolStep = 0.02
MaskRate = 0.025
mkidForceCount = 10
rejectRate = 2
#Shared Memories
bs_start = 0.0
bs_end = 0.0
lockin_start = 0.0
lockin_end = 0.0

def linear_fmin(func,start,end,step):
    min = func(start)
    xmin = start
    points = numpy.arange(start,end,step)
    lenth = len(points)
    div = lenth*1.0 / 100
    nowProg = 0
    print("Start seeking minimum value of function by linear sweeping..")
    print("Start = {0:.4e}, End = {1:.4e}, Step = {2:.4e}".format(start,end,step))
    with progbar.progbar(lenth) as pb:
        for i in range(pb.period()):
            val = func(points[i])
            if min > val:
                min = val
                xmin = points[i]
            prog = math.floor(i/div)
            if prog > nowProg:
                nowProg = prog;
                pb.update(i)
    print("Result = ", xmin)
    return xmin

def findSE(arr2d,freq,timeSpanRate):
    start = -1.0;
    for i in range(0,len(arr2d)):
        if(arr2d[i][1] == arr2d[i+1][1]) and (-arr2d[i][0] + arr2d[i+1][0]) * freq < timeSpanRate:
            start = arr2d[i][0]
            break
    end = -1.0
    for i in range(0,len(arr2d)):
        i = len(arr2d) - i - 2
        if(arr2d[i][1] == arr2d[i+1][1]) and (-arr2d[i][0] + arr2d[i+1][0]) * freq < timeSpanRate:
            end = arr2d[i+1][0]
            break
    return (start,end)

def isInMask(point, bsArr, bsFreq, MaskRate):
    if len([x for x in bsArr if abs(point - x[0]) * bsFreq < MaskRate]) == 0:
        return False
    return True

def lockinBS(bsFile, oppsFile, force=False):
    global bs_start
    global bs_end
    global lockin_start
    global lockin_end
    bsArr_new = []
    if not force and os.path.exists(bsFile+".calib"):
        print("Calib file detected!")
        bsArr_new = numpy.loadtxt(bsFile + ".calib")
    else:
        bsArr = numpy.loadtxt(bsFile)
        oppsArr = numpy.loadtxt(oppsFile)
    
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
    bs_start,bs_end = findSE(bsArr_new, BSFreq, BSSE_TimeSpanRate)
    if bs_end > bs_start:
        print("Found BS start point and end point")
        print("BS_Start = ", bs_start)
        print("BS_End = ", bs_end)
    else:
        print("Can't find BS start point and end point")
        return
    lockin_start,lockin_end = findSE(bsArr_new, BSFreq, LSE_TimeSpanRate)
    lockin_start = math.ceil(lockin_start / lockin_dt) * lockin_dt
    lockin_end = math.floor(lockin_end / lockin_dt) * lockin_dt
    
    return (bsArr_new,scipy.interpolate.interp1d(bsArr_new[:,0],bsArr_new[:,1]))

def lockin(bsArr,interp_BS, mkidFile, force = False):
    mkidArr = []
    mkidShift = 0.0
    print("Lockin File:",mkidFile,flush = True)
    if not force and os.path.exists(mkidFile + ".calib"):
        print("Calib file detected!",flush = True)
        mkidArr = numpy.loadtxt(mkidFile + ".calib")
        mkidShift = mkidArr[0][0]
    else:
        mkidArr = numpy.loadtxt(mkidFile)
        mkidArr = numpy.vstack((mkidArr[:,0],mkidArr[:,2])).T
        #Get calibrated mkid shift
        mkidArr = mkidArr * [1.0/FPSFreq,1.0]
        #print("start,end step=",bs_start,bs_end,ConvolStep)
        mkidArr_mean0 = mkidArr- [0,numpy.mean(mkidArr[:,1])]
        print("Create interpolate function",flush=True)
        interp_mkid_mean0 = scipy.interpolate.interp1d(mkidArr_mean0[:,0],mkidArr_mean0[:,1])
        print("Done",flush=True)
        #def S(A):
        
            #correlate:
            #def F(x):
            #    return interp_mkid(x-A) * interp_BS(x)
            #i1 = scipy.integrate.quad(F,bs_start,bs_end)
            #return i1
        #mkidShift = linear_fmin(S,mkidDefaultShift,mkidDefaultShift+mkidShiftErr,mkidShiftStep)
        #print("Calibrating minimum point by N-M method",flush=True)
        #mkidShift_calib = scipy.optimize.fmin(S,mkidShift)
        #print("Mkid shift time is:",mkidShift_calib)

        #New correlate method:
        bs_start2 = max(mkidArr[0][0]+mkidShiftErr+mkidDefaultShift,bs_start)
        bs_end2 = min(mkidArr[-1][0]+mkidDefaultShift,bs_end)

        xpoints_mkid = numpy.arange(bs_start2-mkidShiftErr-mkidDefaultShift,\
            bs_end2-mkidDefaultShift,mkidShiftStep)
    
        xpoints_BS = numpy.arange(bs_start2,bs_end2,mkidShiftStep)
        corr1 = list(map(interp_mkid_mean0,xpoints_mkid))
        corr2 = list(map(lambda x: -interp_BS(x) + 0.5,xpoints_BS))
        print("Begin correlating, points_mkid = {0:d}, points_bs = {1:d}".format(len(xpoints_mkid),len(xpoints_BS)),flush = True)

        corrArr = scipy.signal.correlate(corr2,corr1,mode = 'valid')
        with open(mkid_file+".corr","wb") as ofs:
            numpy.savetxt(ofs,numpy.vstack((list(map(lambda i:mkidDefaultShift+mkidShiftStep*i,range(0,len(corrArr)))),corrArr)).T)
        mkidShift = corrArr.argmax()* mkidShiftStep + mkidDefaultShift
     
        mkidArr = mkidArr + [mkidShift,0]
        with open(mkidFile + ".calib","wb") as ofs:
            numpy.savetxt(ofs,mkidArr)
        print("Output calibrated mkid data done!")
    print("Mkid shift = ",mkidShift,flush = True)
    #Lock in
    if not force and os.path.exists(mkidFile+".lockin"):
        print("Lockin file detected, skip",flush = True)
    else:
        lockin_start2 = math.ceil(max(lockin_start,mkidArr[0][0]))
        lockin_end2 = math.floor(min(lockin_end,mkidArr[-1][0]))
        max_iter2 = math.floor((lockin_end - lockin_start) / lockin_dt + 0.5)

        iter = 1
        sum_ON = 0.0
        num_ON = 0
        sum_OFF = 0.0
        num_OFF = 0.0
        lockin_Arr = []
        mkidFiltered = [x for x in mkidArr if x[0] >= lockin_start2]
        print("Begin lockin:",flush=True)
        with progbar.progbar(max_iter2) as pb:
            for point in mkidFiltered:
                if point[0] < lockin_start2 + iter * lockin_dt:
                    if isInMask(point[0], bsArr, BSFreq, MaskRate):
                        continue
                    if interp_BS(point[0]) > 0.5:
                        #ON
                        sum_ON = sum_ON + point[1]
                        num_ON = num_ON + 1
                    else:
                        #OFF
                        sum_OFF = sum_OFF + point[1]
                        num_OFF = num_OFF + 1
                else:
                    if num_ON > 0 and num_OFF > 0:
                        lockin_Arr.append([lockin_start2 + (iter - 0.5) * lockin_dt, - sum_ON/num_ON + sum_OFF/num_OFF])
                        sum_ON,sum_OFF = 0.0, 0.0
                        num_ON,num_OFF = 0, 0
                        iter = iter + 1
                        pb.update(iter)
                        if(iter > max_iter2):
                            break
                    else:
                        print("Illigal data has occured!")
                        return
        with open(mkidFile+".lockin","wb") as ofs:
            numpy.savetxt(ofs,numpy.array(lockin_Arr))
        print("Output lockin data done!",flush=True)
    return mkidShift

if not len(sys.argv) < 4:
    force = False
    mkid_file_count = len(sys.argv[3:])
    if mkid_file_count <= mkidForceCount:
        force = True
    bsArr,interp_BS = lockinBS(sys.argv[1],sys.argv[2],force)
    shiftlist = []
    startTime = time.time()
    
    for mkid_file in sys.argv[3:]:
        mkidNoRe = re.compile("MKID(?P<No>\d\d\d)")
        mkidNoMatch = mkidNoRe.search(mkid_file)
        mkidNo = int(mkidNoMatch.group('No'))
        print("MKID No:",mkidNo,flush=True)
        shift = lockin(bsArr,interp_BS,mkid_file,force)
        print("Shift = ",shift,flush = True)
        shiftlist.append([float(mkidNo),shift])
    shiftArr = numpy.array(shiftlist)
    #shiftArr.sort(0)
    if not force:
        with open("shiftlist.txt","wb") as ofs:
            numpy.savetxt(ofs,shiftArr[shiftArr[:,0].argsort(),:])
        print("Shift list has generated!")
        m = numpy.mean(shiftArr[:,1])
        s = numpy.std(shiftArr[:,1])
        shiftArr_rejected = [x for x in shiftArr if abs(x[1] - m) >= s*rejectRate]
        if len(shiftArr_rejected) > 0:
            with open("rejectlist.txt","wb") as ofs:
                numpy.savetxt(ofs,shiftArr_rejected)
            print("{0:d} points are illegal!".format(len(shiftArr_rejected)))
        print("Total used time is {0:d} seconds".format(int(time.time()-startTime)))
