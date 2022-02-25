import numpy as nmp
def GandBbus(BusNo,Ybus):
    Vbase = 230.94
    Sbase = 10e3
    Zbase = Vbase**2/Sbase
    nphase = 4
    NoColY = BusNo*nphase

    NoElem = len(Ybus)
    GRowCount =0
    GColCount =0
    BRowCount =0
    BColCount =0
    B=[ [0]*NoColY for i in range(NoColY)]
    G=[ [0]*NoColY for i in range(NoColY)]
    for NoElem in range(NoElem):
        if NoElem%2==0:
            if abs(Ybus[NoElem])>100000:
                G[GRowCount][GColCount]=2666.67*Zbase
            else:
                G[GRowCount][GColCount]=Ybus[NoElem]*Zbase
            GColCount +=1
            if GColCount ==NoColY:
                GColCount=0
                GRowCount +=1
        else:
            if abs(Ybus[NoElem])>100000:
                B[BRowCount][BColCount]=-1166.666667*Zbase
            else:
                B[BRowCount][BColCount]=Ybus[NoElem]*Zbase
            BColCount +=1
            if BColCount ==NoColY:
                BColCount=0
                BRowCount +=1
    return G,B

def Index_calc(Bus1,Bus2,Ph1,Ph2):
    start1 = Bus1*4
    num1 = start1+Ph1
    start2 = Bus2*4
    num2 = start2+Ph2
    return num1,num2

def Index_calcJ(Bus1,Bus2,Ph1,Ph2):
    start1= Ph1*75
    num1 = start1+Bus1
    start2 = Bus2*4
    num2 = start2+Ph2
    return num1,num2

def Index_calcJ2(Bus,Ph):
    start = Bus*4
    num = start+Ph
    return num

def iter_reduction(G,B):
    non_zero = nmp.where(nmp.absolute(G+1j*B)>0)
    iter_nonzero = [[0]*4 for i in range(len(non_zero[0]))]
    for i in range(len(non_zero[0])):
        i = int(i)
        ph_detect1 = non_zero[0][i]%4
        ph_detect2 = non_zero[1][i]%4
        Bus1 = (non_zero[0][i]-ph_detect1)/4
        Bus2 = (non_zero[1][i]-ph_detect2)/4
        iter_nonzero[i] = [int(Bus1),int(Bus2),int(ph_detect1),int(ph_detect2)]
    return iter_nonzero