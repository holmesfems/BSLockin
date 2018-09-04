#!/usr/bin/env python3
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
import param
#Shared Memories
bs_start = 0.0
bs_end = 0.0
lockin_start = 0.0
lockin_end = 0.0
Param = param.Param.copy()

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

def isInMask(point, bsArr, bsfreq, maskrate):
    if (numpy.abs(bsArr - point) < (maskrate / bsfreq)).any():
        return True
    return False

def lockinBS(bsFile, oppsFile):
    global bs_start
    global bs_end
    global lockin_start
    global lockin_end
    bsArr_new = []
    if not Param['force'] and os.path.exists(bsFile+".calib"):
        print("Calib file detected!")
        bsArr_new = numpy.loadtxt(bsFile + ".calib")
    else:
        bsArr = numpy.loadtxt(bsFile)
        oppsArr = numpy.loadtxt(oppsFile)
    
        #Calculate calibrated DAQ Sampling Frequency
        firstPulse = -1
        for i in range(0,len(oppsArr)):
            if (oppsArr[i][1] == 1) and (oppsArr[i+1][1] == 1) and (( -oppsArr[i][0] + oppsArr[i+1][0]) / Param['DAQFreq'] < 0.1):
                firstPulse = oppsArr[i][0]
                break
        lastPulse = -1
        for i in range(0,len(oppsArr)):
            i = len(oppsArr) - i - 2
            if (oppsArr[i][1] == 1) and (oppsArr[i+1][1] == 1) and (( -oppsArr[i][0] + oppsArr[i+1][0]) / Param['DAQFreq'] < 0.1):
                lastPulse = oppsArr[i][0]
                break
        if firstPulse < 0 or lastPulse < 0:
            print("Failed to fit opps pulse!")
            return
        span = round((lastPulse-firstPulse)/Param['DAQFreq'])
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
    bs_start,bs_end = findSE(bsArr_new, Param['BSFreq'], Param['BSSE_TimeSpanRate'])
    if bs_end > bs_start:
        print("Found BS start point and end point")
        print("BS_Start = ", bs_start)
        print("BS_End = ", bs_end)
    else:
        print("Can't find BS start point and end point")
        return
    lockin_start,lockin_end = findSE(bsArr_new, Param['BSFreq'], Param['LSE_TimeSpanRate'])
    lockin_start = math.ceil(lockin_start / Param['lockin_dt']) * Param['lockin_dt']
    lockin_end = math.floor(lockin_end / Param['lockin_dt']) * Param['lockin_dt']
    
    return (bsArr_new,scipy.interpolate.interp1d(bsArr_new[:,0],bsArr_new[:,1]))

def lockin(bsArr,interp_BS, mkidFile, shiftTime = None):
    mkidArr = []
    mkidShift = 0.0
    print("Lockin File:",mkidFile,flush = True)
    doCalib = True
    corr_isPositive = None
    #generate shift time
    if not Param['force'] and os.path.exists(mkidFile + ".calib"):
        print("Calib file detected!",flush = True)
        mkidArr = numpy.loadtxt(mkidFile + ".calib")
        mkidShift = mkidArr[0][0]
        if mkidShift == shiftTime and os.path.exists(mkidFile + ".corr"):
            doCalib = False
            corrArr = numpy.loadtxt(mkidFile + ".corr")
            corr_isPositive = (scipy.interpolate.interp1d(corrArr[:,0],corrArr[:,1])(mkidShift) > 0)
        else:
            print("Use detected shift time:",shiftTime,flush = True)

    elif not shiftTime == None:
        print("Use detected shift time:",shiftTime,flush = True)
        mkidShift = shiftTime

    
    if doCalib:
        mkidArr = numpy.loadtxt(mkidFile)[:,0:3:2] / [Param['FSPFreq'],1]
        #Get calibrated mkid shift
        #print("start,end step=",bs_start,bs_end,Param['ConvolStep'])
        mean = numpy.mean(mkidArr[:,1])
        mkidArr_mean0 = mkidArr - [0,mean]
        
        print("Create interpolate function",flush=True)
        interp_mkid_mean0 = scipy.interpolate.interp1d(mkidArr_mean0[:,0],mkidArr_mean0[:,1],fill_value = (mkidArr_mean0[0][1],mkidArr_mean0[-1][1]))
        print("Done",flush=True)
        _mkidDefaultShift = 0.0
        #New correlate method:
        if not shiftTime == None:
            _mkidShiftErr = Param['mkidShiftErr_useShift']
            _mkidDefaultShift = shiftTime - 0.5 * _mkidShiftErr + Param['mkidBaseShift']
        else:
            _mkidShiftErr = Param['mkidShiftErr']
            _mkidDefaultShift = Param['mkidDefaultShift']
        bs_start2 = 0.0
        bs_end2 = 0.0
        shiftstart = 0.0
        shiftend = 0.0
        if Param['use_costumSE']:
            print("use costum SE: %.3f ~ %.3f" % (Param['mkidShift_costumS'], Param['mkidShift_costumE']))
            shiftstart = (1 - Param['mkidShift_costumS']) * mkidArr_mean0[0][0] + (Param['mkidShift_costumS']) * mkidArr_mean0[-1][0] + _mkidShiftErr
            shiftend = (1 - Param['mkidShift_costumE']) * mkidArr_mean0[0][0] + (Param['mkidShift_costumE']) * mkidArr_mean0[-1][0] + _mkidShiftErr
        else:
            shiftstart = bs_start
            shiftend = bs_end
        bs_start2 = max(mkidArr_mean0[0][0]+_mkidShiftErr+_mkidDefaultShift,shiftstart)
        bs_end2 = min(mkidArr_mean0[-1][0]+_mkidDefaultShift ,shiftend)

        xpoints_mkid = numpy.arange(bs_start2-_mkidShiftErr-_mkidDefaultShift,\
            bs_end2-_mkidDefaultShift,Param['mkidShiftStep'])
    
        xpoints_BS = numpy.arange(bs_start2,bs_end2,Param['mkidShiftStep'])
        corr1 = interp_mkid_mean0(xpoints_mkid)
        corr2 = -interp_BS(xpoints_BS)+0.5

        print("Begin correlating, points_mkid = {0:d}, points_bs = {1:d}".format(len(xpoints_mkid),len(xpoints_BS)),flush = True)

        corrArr = scipy.signal.correlate(corr2,corr1,mode = 'valid')
        with open(mkid_file+".corr","wb") as ofs:
            numpy.savetxt(ofs,numpy.array([[_mkidDefaultShift+x*Param['mkidShiftStep'],corrArr[x]] for x in range(0,len(corrArr))]))
        argmax = numpy.abs(corrArr).argmax()
        mkidShift = argmax * Param['mkidShiftStep'] + _mkidDefaultShift
        if corrArr[argmax] < 0:
            corr_isPositive = False
            print("Warning: detected minus correlation function")
        else:
            corr_isPositive = True
        mkidArr = mkidArr_mean0 + [mkidShift,0]
        with open(mkidFile + ".calib","wb") as ofs:
            numpy.savetxt(ofs,mkidArr+[0,mean])
        print("Output calibrated mkid data done!")
        Param['force'] = True
    print("Mkid shift = ",mkidShift,flush = True)
    #Lock in
    if not Param['force'] and os.path.exists(mkidFile+".lockin"):
        print("Lockin file detected, skip",flush = True)
    else:
        lockin_start2 = math.ceil(max(lockin_start,mkidArr[0][0]))
        lockin_end2 = math.floor(min(lockin_end,mkidArr[-1][0]))
        max_iter2 = math.floor((lockin_end2 - lockin_start2) / Param['lockin_dt'] + 0.5)

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
                if not point[0] < lockin_start2 + iter * Param['lockin_dt']:
                    if num_ON > 0 and num_OFF > 0:
                        lockin_Arr.append([lockin_start2 + (iter - 0.5) * Param['lockin_dt'], - sum_ON/num_ON + sum_OFF/num_OFF if corr_isPositive == None or corr_isPositive else sum_ON/num_ON - sum_OFF/num_OFF])
                        sum_ON,sum_OFF = 0.0, 0.0
                        num_ON,num_OFF = 0, 0
                        iter = iter + 1
                        pb.update(iter)
                        if(iter >= max_iter2):
                            break
                    else:
                        print("Illigal data has occured!")
                        return
                if isInMask(point[0], bsArr[:,0], Param['BSFreq'], Param['MaskRate']):
                    continue
                if interp_BS(point[0]) > 0.5:
                    #ON
                    sum_ON = sum_ON + point[1]
                    num_ON = num_ON + 1
                else:
                    #OFF
                    sum_OFF = sum_OFF + point[1]
                    num_OFF = num_OFF + 1
        with open(mkidFile+".lockin","wb") as ofs:
            numpy.savetxt(ofs,numpy.array(lockin_Arr))
        print("Output lockin data done!",flush=True)
    return (mkidShift,corr_isPositive)

paramSetRe = re.compile(r"(?P<param>[^=]+)=(?P<value>[^=]+)")
paramSetList = [x for x in sys.argv if not paramSetRe.match(x) is None]
intRe = re.compile(r"-?\d+")
boolRe = re.compile(r"(True)|(False)")
#set params
for item in paramSetList:
    match = paramSetRe.match(item)
    paramName=match.group('param')
    if not paramName in Param:
        continue
    value=match.group('value')
    if not intRe.match(value) is None:
        Param[paramName] = int(value)
    elif not boolRe.match(value) is None:
        Param[paramName] = True if value == 'True' else False
    else:
        Param[paramName] = float(value)
    print("Set param %s to %s"%(paramName,value))
print("now Param:")
print(Param)

FilteredParam = [x for x in sys.argv if not x in paramSetList]
if not len(FilteredParam) < 4:

    mkid_file_count = len(FilteredParam[3:])
    if mkid_file_count <= Param['mkidForceCount']:
        Param['force'] = True
    bsArr,interp_BS = lockinBS(FilteredParam[1],FilteredParam[2])
    shiftlist = []
    shiftReadlist = {}
    corrNegativelist = []
    startTime = time.time()
    if Param['use_shift']:
        if(os.path.exists("shiftlist.txt")):
            shiftReadlist = {math.floor(x[0]+0.5):x[1] for x in numpy.loadtxt("shiftlist.txt")}
            print("Detected shiftlist",flush=True)
        else:
            Param['use_shift'] = False
    for mkid_file in FilteredParam[3:]:
        mkidNoRe = re.compile("MKID(?P<No>\d\d\d)")
        mkidNoMatch = mkidNoRe.search(mkid_file)
        mkidNo = int(mkidNoMatch.group('No'))
        print("MKID No:",mkidNo,flush=True)
        shift,corr_isPositive = None,None
        if Param['use_shift']:
            shift,corr_isPositive = lockin(bsArr,interp_BS,mkid_file,shiftTime = shiftReadlist[mkidNo])
        else:
            shift,corr_isPositive = lockin(bsArr,interp_BS,mkid_file)
        #print("Shift = ",shift,flush = True)
        if not shift == None:
            shiftlist.append([float(mkidNo),shift])
        if not corr_isPositive:
            corrNegativelist.append(mkidNo)
    shiftArr = numpy.array(shiftlist)
    #shiftArr.sort(0)
    shiftlistName = "shiftlist2.txt" if Param['use_shift'] else "shiftlist.txt"
    if not mkid_file_count <= Param['mkidForceCount']:
        if len(shiftArr) > 0:
            with open(shiftlistName,"wb") as ofs:
                numpy.savetxt(ofs,shiftArr[shiftArr[:,0].argsort(),:])
            print("Shift list has generated!")
            m = numpy.mean(shiftArr[:,1])
            s = numpy.std(shiftArr[:,1])
            shiftArr_rejected = [x for x in shiftArr if abs(x[1] - m) >= s*Param['rejectRate']]
            if len(shiftArr_rejected) > 0:
                with open("rejectlist.txt","wb") as ofs:
                    numpy.savetxt(ofs,shiftArr_rejected)
                print("{0:d} points are illegal!".format(len(shiftArr_rejected)))
                strrejectx = ["%03d" % math.floor(x[0]+0.5) for x in shiftArr_rejected]
                print("Illigal mkid list: "+", ".join(strrejectx))
        if len(corrNegativelist) > 0:
            print("{0:d} mkids has negative correlation function!".format(len(corrNegativelist)))
            corrNegativelistStr = ["%03d" % x for x in corrNegativelist]
            print("Negative mkid list: ",", ".join(corrNegativelistStr))
            with open("negativelist.txt","w") as ofs:
                ofs.write("\n".join(corrNegativelistStr))
    print("Total used time is {0:d} seconds".format(int(time.time()-startTime)))

else:
    print("Usage: lockin.py [BS_DATA.out] [OPPS_DATA.out] [MKID_DATA]+")
