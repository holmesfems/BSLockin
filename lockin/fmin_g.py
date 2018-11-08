import math
def fmin_g(func,start,end,acc,L=None,R=None,ML=None,MR=None):
    """
    Search fmin location by golden-section search
    """
    if abs(start-end) < abs(acc):
        xM = 0.5*start+0.5*end
        return (xM,start,end)
    phi = (math.sqrt(5.0)-1.0)/2.0
    xML = phi*start+(1-phi)*end
    xMR = phi*end+(1-phi)*start
    if (bool(L is None)!=bool(R is None)):
        #Error occured
        #print("error occured")
        return(None,start,end)
    if L is None:
        L = func(start)
    if R is None:
        R = func(end)
    if ML is None:
        ML = func(xML)
    def fmin_L():
        return fmin_g(func,start,xMR,acc,L=L,R=MR,MR=ML)
    def fmin_R():
        return fmin_g(func,xML,end,acc,L=ML,R=R,ML=MR)
    if L<=ML and R <=ML:
        xM = 0.5*start+0.5*end
        print("L={0} R={1} ML={2}".format(L,R,ML))
        return (xM,start,end)
    if MR is None:
        MR = func(xMR)
    if MR > ML:
        return fmin_L()
    else:
        return fmin_R()

if __name__ == '__main__':
    print(fmin_g(lambda x:1/(x+1.1)+x,-1,0.5,0.01))