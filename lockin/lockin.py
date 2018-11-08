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
import fmin_g
#Shared Memories
bs_start = 0.0
bs_end = 0.0
lockin_start = 0.0
lockin_end = 0.0
Param = param.Param.copy()

def findSE(arr2d,freq,timeSpanRate):
    diffX = numpy.diff(arr2d[:,0])
    diffY = numpy.diff(arr2d[:,1])
    seLogic = numpy.logical_and(diffY == 0,diffX*freq < timeSpanRate)
    start = arr2d[:-1][seLogic][0][0]
    end = arr2d[1:][seLogic][-1][0]
    return (start,end)

def rms(arr1d):
    mean=numpy.mean(arr1d)
    rms=math.sqrt(numpy.sum((arr1d-mean)**2.0)/len(arr1d))
    return rms

def genSplitIndex(arr,area):
    return numpy.bincount(numpy.digitize(arr,area)).cumsum()[:-1]

def genMaskedBS(bsArr):
    maskrate = Param['MaskRate']
    bsfreq = Param['BSFreq']
    quant = 1.0/Param['DAQFreq']
    mrdt = maskrate/bsfreq
    bsm = bsArr[1:-1][numpy.diff(bsArr[1:,1])!=0]
    bsp = bsArr[1:-1][numpy.diff(bsArr[1:,1])==0]
    bsmm = bsm + [-mrdt-0.5*quant,0]
    bsmp = (bsm + [-mrdt+0.5*quant,0])*[1,0]+[0,0.5]
    bspm = (bsp + [+mrdt-0.5*quant,0])*[1,0]+[0,0.5]
    bspp = bsp + [+mrdt+0.5*quant,0]
    res = numpy.vstack((bsArr[0],numpy.hstack((bsmm,bsmp,bspm,bspp)).reshape(-1,2),bsArr[-1]))
    return res

def genMaskedMkid(mkidArr,bsInterp):
    return mkidArr[numpy.abs((bsInterp(mkidArr[:,0])-0.5))>=0.4]

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
        lastPulse = -1
        diffX = numpy.diff(oppsArr[:,0])
        diffY = numpy.diff(oppsArr[:,1])
        seLogic = numpy.logical_and(oppsArr[:-1,1]==1,diffY==0,diffX/Param['DAQFreq']<0.1)
        flPulse = oppsArr[:-1][seLogic]
        if len(flPulse) >= 2:
            firstPulse = flPulse[0][0]
            lastPulse = flPulse[-1][0]
        if firstPulse < 0 or lastPulse < 0:
            print("Failed to fit opps pulse!")
            return
        span = round((lastPulse-firstPulse)/Param['DAQFreq'])
        DAQFreq_calib = (lastPulse - firstPulse) / span
    
        #Calculate Zero point of time:
        #Zero Point is firstPulse of opps
        #Get first timestamp of BS
        bsFileOpen = open(bsFile)
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
    return bsArr_new

def lockin(bsArr,interp_BS, mkidFile, shiftTime = None, bsArr_masked = None, interp_BS_masked = None):
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
        #bs_start2 = 0.0
        #bs_end2 = 0.0
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
        corr2 = -interp_BS_masked(xpoints_BS)+0.5 if Param['shift_withMask'] else -interp_BS(xpoints_BS)+0.5

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
    lockin_Arr = []
    if not Param['force'] and os.path.exists(mkidFile+".lockin"):
        print("Lockin file detected, skip",flush = True)
    else:
        lockin_start2 = math.ceil(max(lockin_start,mkidArr[0][0]))
        lockin_end2 = math.floor(min(lockin_end,mkidArr[-1][0]))
        max_iter2 = math.floor((lockin_end2 - lockin_start2) / Param['lockin_dt'] + 0.5)
        mkidFiltered = mkidArr[mkidArr[:,0]>lockin_start2]
        mkidOutofMask = genMaskedMkid(mkidFiltered,interp_BS_masked)
        print("Begin lockin:",flush=True)
        with progbar.progbar(max_iter2) as pb:
            lockinArea = numpy.arange(lockin_start2+Param['lockin_dt'],lockin_end2,Param['lockin_dt'])
            mkidOutofMask_split = numpy.split(mkidOutofMask,genSplitIndex(mkidOutofMask[:,0],lockinArea))[:-1]
            split_ON = [x[interp_BS(x[:,0])>0.5] for x in mkidOutofMask_split]
            pb.update(max_iter2/4)
            split_OFF = [x[interp_BS(x[:,0])<=0.5] for x in mkidOutofMask_split]
            pb.update(max_iter2/2)
            sum_ON = numpy.array([numpy.sum(x[:,1]) for x in split_ON])
            sum_OFF = numpy.array([numpy.sum(x[:,1]) for x in split_OFF])
            count_ON = numpy.array([len(x) for x in split_ON])
            count_OFF = numpy.array([len(x) for x in split_OFF])
            res = sum_ON / count_ON - sum_OFF / count_OFF
            if len(lockinArea)!=len(res):
                lockinArea = lockinArea[:len(res)]
            lockin_Arr = numpy.vstack((lockinArea-0.5*Param['lockin_dt'],res)).T
        with open(mkidFile+".lockin","wb") as ofs:
            numpy.savetxt(ofs,numpy.array(lockin_Arr))
        with open(mkidFile+".outmask","wb") as ofs:
            numpy.savetxt(ofs,numpy.array(mkidOutofMask))
        print("Output lockin data done!",flush=True)
    return (mkidShift,corr_isPositive,lockin_Arr)

#Main process
paramSetRe = re.compile(r"(?P<param>[^=]+)=(?P<value>[^=]+)")
paramSetList = [x for x in sys.argv if not paramSetRe.match(x) is None]
intRe = re.compile(r"-?\d+$")
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
    bsArr = lockinBS(FilteredParam[1],FilteredParam[2])
    interp_BS = scipy.interpolate.interp1d(bsArr[:,0],bsArr[:,1])
    shiftlist = []
    shiftReadlist = {}
    corrNegativelist = []
    mrList = []
    startTime = time.time()
    if Param['use_shift']:
        if(os.path.exists("shiftlist.txt")):
            shiftReadlist = {math.floor(x[0]+0.5):x[1] for x in numpy.loadtxt("shiftlist.txt")}
            print("Detected shiftlist",flush=True)
        else:
            Param['use_shift'] = False
    for mkid_file in FilteredParam[3:]:
        mkidNoRe = re.compile(r"MKID(?P<No>\d\d\d)")
        mkidNoMatch = mkidNoRe.search(mkid_file)
        mkidNo = int(mkidNoMatch.group('No'))
        print("MKID No:",mkidNo,flush=True)
        shift,corr_isPositive = None,None
        if Param['fit_mr']:
            def rms_mask(MR):
                global shift,corr_isPositive
                Param['MaskRate'] = MR
                bsArr_masked = genMaskedBS(bsArr)
                interp_BS_masked = scipy.interpolate.interp1d(bsArr_masked[:,0],bsArr_masked[:,1])
                shift,corr_isPositive,lockin_Arr = lockin(bsArr,interp_BS,mkid_file,\
                    shiftTime = shiftReadlist[mkidNo] if Param['use_shift'] else None,\
                    bsArr_masked = bsArr_masked,\
                    interp_BS_masked = interp_BS_masked\
                    )
                rms_s = math.ceil(len(lockin_Arr)*Param['rms_s'])
                rms_e = math.floor(len(lockin_Arr)*Param['rms_e'])
                rms_val = rms(lockin_Arr[rms_s:rms_e])
                return rms_val
            mr,st,ed = fmin_g.fmin_g(rms_mask,Param['mr_s'],Param['mr_e'],Param['mr_a'])
            Param['MaskRate'] = mr
            print("Masking Rate fitting result:{0}".format(mr))
            print("st={0} ed={1} ac={2}".format(st,ed,ed-st))
            if ed-st > Param['mr_a']:
                print("Warning: ac is larger than expected")
            mrList.append([float(mkidNo),mr,ed-st])
        bsArr_masked = genMaskedBS(bsArr)
        interp_BS_masked = scipy.interpolate.interp1d(bsArr_masked[:,0],bsArr_masked[:,1])
        shift,corr_isPositive,lockin_Arr = lockin(bsArr,interp_BS,mkid_file,\
            shiftTime = shiftReadlist[mkidNo] if Param['use_shift'] else None,\
            bsArr_masked = bsArr_masked,\
            interp_BS_masked = interp_BS_masked\
            )
        #print("Shift = ",shift,flush = True)
        if not shift == None:
            shiftlist.append([float(mkidNo),shift])
        if not corr_isPositive:
            corrNegativelist.append(mkidNo)
    mrList = numpy.array(mrList)
    shiftArr = numpy.array(shiftlist)
    #shiftArr.sort(0)
    shiftlistName = "shiftlist2.txt" if Param['use_shift'] else "shiftlist.txt"
    if not mkid_file_count <= Param['mkidForceCount']:
        if len(mrList) > 0:
            with open("mrlist.txt","wb") as ofs:
                numpy.savetxt(ofs,mrList[mrList[:,0].argsort()])
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
