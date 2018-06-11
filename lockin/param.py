#Params
lockin_dt = 0.25
mkidDefaultShift = 2.0
BSFreq = 10.0
DAQFreq = 80000.0
FSPFreq = 128.0
BSSE_TimeSpanRate = 10
LSE_TimeSpanRate = 0.55
mkidBaseShift = -0.02
mkidShiftErr = 4.0
mkidShiftErr_useShift = 0.05
mkidShift_costumS = 1.0*(25 - 15)/(80 - 15)
mkidShift_costumE = 1.0*(35 - 15)/(80 - 15)
mkidShiftStep = 0.001
#mkidShiftStep_useShift = 0.0005
ConvolStep = 0.02
MaskRate = 0.025
mkidForceCount = 10
rejectRate = 2

force = True
use_shift = True
use_costumSE = True