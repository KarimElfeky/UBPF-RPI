def Load_extract(Busnames,Loadbuses,Loadnames,Loadshapenames,dssLoad,dssLoadshapes):
    BusNo = len(Busnames)
    P = [ [0]*3 for i in range(BusNo)]
    Q = [ [0]*3 for i in range(BusNo)]
    NoLoads = len(Loadnames)
    Buscount = 0
    loadphase = 0
    for Loadcount in range(NoLoads):
        dssLoad.Name = Loadnames[Loadcount]
        if Loadnames[Loadcount] in Loadshapenames:
            dssLoadshapes.Name = Loadnames[Loadcount]
            LB = Loadbuses[Buscount]
            busind = Busnames.index(str(LB))
            P[busind][loadphase]=float(dssLoadshapes.Pmult[0])*float(dssLoad.kW)
            Q[busind][loadphase]=float(dssLoadshapes.Qmult[0])*float(dssLoad.kvar)
            loadphase+=1
            if loadphase==3:
                Buscount +=1
                loadphase=0
        else:
            LB = str(76)
            busind = Busnames.index(str(LB))
            P[busind][loadphase]=P[busind][loadphase]+float(dssLoad.kW)
            Q[busind][loadphase]=Q[busind][loadphase]+float(dssLoad.kvar)
            loadphase +=1
            if loadphase==3:
                Buscount +=1
                loadphase=0
    return P,Q
