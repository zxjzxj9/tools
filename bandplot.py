#! /usr/bin/python 

from lxml import etree
import matplotlib.pyplot as plt
import numpy as np
from math import floor,acos
import itertools
from scipy.interpolate import interp1d,Rbf,\
                              InterpolatedUnivariateSpline
from scipy import interp
import sys
#from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

tree = etree.parse("vasprun.xml")
root = tree.getroot()

kpname = ["$\Gamma$","F","Q","Z","B","$\Gamma$"]

Emin = -5.0 # eV
Emax =  5.0  # eV

interp = True
dosplot = False

# init data
upper = 5 # upper limit of VB
lower = 5 # lower limit of CB

kplist = []
kplistdiv = []
kpdiv = 0
lattice = []
reclattice = []
efermi = 0.0
spin1 = {}
spin2 = {}

dos1 = []
dos2 = []

def parseKpoint(pn):
    data = []
    div = 0
    for line in pn:
        if line.tag == 'i':
            div = int(line.text)
        elif line.tag == 'v':
            data.append(map(lambda x: float(x), line.text.split()))
    return (div,data)

def parseVal(pn,ldos=False):
    data = []
    for line in pn:
        if line.tag == 'v':
            data.append(map(lambda x: float(x), line.text.split()))
        if line.tag == 'r' and not ldos:
            data.append(float(line.text.split()[0]))
        if line.tag == 'r' and ldos:
            data.append((float(line.text.split()[0]),\
                        float(line.text.split()[1])))
    return data

def parseBand(pn):
    spin1 = {}
    spin2 = {}
    for setdata in pn:
        if setdata.get("comment") == "spin 1":
            for snode in setdata:
                number = int(snode.get("comment").split()[1])
                spin1[number] = parseVal(snode)
        if setdata.get("comment") == "spin 2":
            for snode in setdata:
                number = int(snode.get("comment").split()[1])
                spin2[number] = parseVal(snode)
    return (spin1,spin2)

for node in root:
    #print node.tag
    if node.tag == "kpoints":
        for kpnode in node:
            #print dir(kpnode)
            #print kpnode.get("param")
            if kpnode.tag == "generation":
                kpdiv,kplist = parseKpoint(kpnode)
            elif kpnode.tag == "varray" and kpnode.get("name") == "kpointlist":
                kplistdiv = parseVal(kpnode)

    if node.tag == "structure":
        for snode in node:
            if snode.tag == "crystal":
                for ssnode in snode:
                    if ssnode.get("name") == "basis":
                        lattice = parseVal(ssnode)
                    if ssnode.get("name") == "rec_basis":
                        reclattice = parseVal(ssnode)
#                        print "reclattice", reclattice

    if node.tag == "calculation":
        for snode in node:
            if snode.tag == "eigenvalues":
                kpdataline = snode.find(".//set")
                #print snode.getchildren()
                #kpdataline = etree.SubElement(snode,"set")
                #print kpdataline.tag
                spin1, spin2 = parseBand(kpdataline)
                #print spin2
            if snode.tag == "dos":
                #efermi = float(etree.SubElement(snode,"i").text)
                efermi = float(snode.find(".//i").text)
                #print snode.tag
                print "Fermi level obtained, values:",efermi

                if dosplot:
                    dosset = \
                    snode.find(".//total").find(".//array").find(".//set")

                    for sset in dosset:
                        if sset.get("comment")=="spin 1":
                            dos1 = parseVal(sset,ldos=True)
                        if sset.get("comment")=="spin 2":
                            dos2 = parseVal(sset,ldos=True)

newBandMap = {}
newBandMap = spin1
#for n in range(1,len(kplistdiv)+1):
#    bandValList = spin1[n]
#    upperBand = filter(lambda x: x>efermi,bandValList)
#    lowerBand = filter(lambda x: x<efermi,bandValList)
#    upperBand.sort()
#    lowerBand.sort()
#    newBandMap[n] = map(lambda x: x - \
#                    efermi,lowerBand[len(lowerBand)-lower:] + lowerBand[0:upper])
#    #print newBandMap[n]

X = []
#distlist = []
kmatrix = np.array(reclattice)
weight = []
#print kplist
for pair in zip(kplist,kplist[1:]):
#    print pair
    dk = np.array(pair[1]) - np.array(pair[0])
#    print dk
#    print kmatrix
    dist = np.linalg.norm(np.dot(kmatrix, dk))
#   print dist
#    distlist.append(dist)
    weight.append(dist)

dx = [0.0]
for i in range(0,len(kplistdiv)):
    dx.append(weight[int(floor(float(i)/kpdiv))])

dx[0] = 0.0
X.append(0.0)
cnt = 0
for i in range(1,len(kplistdiv)):
#    if i%kpdiv == 0: continue
    X.append(X[cnt] + dx[i])
    cnt += 1
#print len(X)
fig = plt.figure()
plt.xlim(0,max(X))
plt.ylim(Emin,Emax)

labelnum = [0.0]
for i in range(1,len(kpname)):
    labelnum.append(labelnum[i-1] + kpdiv*weight[i-1])


#axeslist = fig.add_axes(labelnum)
#print labelnum

plt.xticks(labelnum,kpname)
plt.ylabel("Energy / eV")
for t in labelnum:
    #print t
    plt.plot([t,t],[Emin,Emax],"k-.")

#for i in range(0,upper+lower):
#    y = []
#    for j in range(1,len(kplistdiv)+1):
#
#    #    if j%kpdiv == 0 and j != len(kplistdiv): continue
#        y.append(newBandMap[j][i])
#    #print len(y)
#    plt.plot(x,y,"b-")
# above code is wrong

def project(n,labeln):
    if n%kpdiv ==0 or n%kpdiv==1:
        return labeln[n/kpdiv] # python 2.0, in 3.0 should be //
    else:
        return labeln[n/kpdiv] + \
        float(n%kpdiv)*(labeln[n/kpdiv+1]-labeln[n/kpdiv])/kpdiv
#    pos = np.array(kplistdiv[n-1])
#    cnt = 0
#    for pair in zip(kplist,kplist[1:]):
#        cnt += 1
#        start = np.array(pair[0])
#        end = np.array(pair[1])
#
#        if within(pos,start,end):
#            return within(pos,start,end) * ( labeln[cnt] - labeln[cnt-1]
#            ) + labeln[cnt-1]

def within(a,b,c):
    # return ratio if a in b c
    if np.linalg.norm(a-b) < 1.0e-3: return 0.0
    if np.linalg.norm(a-c) < 1.0e-3: return 1.0
    if np.linalg.norm(np.cross(a-b,a-c))<1.0e-3: return \
        np.linalg.norm(a-b)/np.linalg.norm(b-c)
    return None

#### Scan the bands next to the Fermi level:
for i in range(len(newBandMap[1])):
    x = []
    y = []
    for j in range(1,len(kplistdiv)+1):
        x.append(project(j,labelnum))
        y.append(newBandMap[j][i] - efermi)
    if interp:
        xnew = []
        ynew = []
        # divide sgements
        cnt = 0
        for pair in zip(labelnum,labelnum[1:]):
            # print pair
            zt = filter(lambda tmp: tmp[0] >= pair[0] and tmp[0] <= pair[1],\
                          zip(x,y))
            sorted(zt, key = lambda x: x[0])
            if zt[0][0] == zt[1][0]: zt.pop(0)
            if zt[-1][0] == zt[-2][0]: zt.pop()
            xt,yt = zip(*zt)

            #print xt
            f = interp1d(xt, yt, kind='cubic')
            #f = InterpolatedUnivariateSpline(xt, yt, k=5)
            #f = Rbf(xt,yt)

            tmpx = np.linspace(pair[0],pair[1],100,endpoint=True)

            xnew.append(tmpx)
            ynew.append(f(tmpx))


            cnt += 1
        x = np.hstack(xnew)
        y = np.hstack(ynew)
        #print xnew
        #print ynew
        #sys.exit()
    #x = np.array(range(len(kplistdiv)))*labelnum[-1]/len(kplistdiv)
    #print x[-1]
    #plt.plot(range(len(x)),x)
    #break
    plt.plot(x,y,"b-")
#    plt.plot(x,y)

plt.show()

if dosplot:
    fig2 = plt.figure()
    #print dos1
    x,y = zip(*dos1)
    x = np.array(x) - efermi
    plt.xlim(Emin,Emax) 
    plt.ylim(0.0,max(zip(*filter(lambda t:t[0]>Emin and t[0]<Emax, \
             dos1))[1])*1.1)
    f = interp1d(x,y,kind='cubic')
    xnew = np.linspace(Emin,Emax,1000)
    ynew = f(xnew)
    #plt.plot(x,y)
    plt.plot(xnew,ynew,color='r')
    plt.show()

