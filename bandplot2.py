"""

Written by zxj, 2016/6/26
Plot the band structure read from vasprun.xml
Version 0.2

Requirements:
lxml, for xml file reading
matplotlib, for plotting
scipy, for interpolation
ConfigParser, to parse ini files

Further extensions:
PySide, add a gui tool 

"""

#! /usr/bin/env python

import os
import sys
import ConfigParser
import argparse
from lxml import etree
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import interp1d
import copy
import scipy as sp
#import lxml
#import pdb


def main_cmd(optFile, spin = 1):
    if not os.path.exists(optFile):
        print("No such configuration file: %s, exiting..." %optFile)
        sys.exit()
    else:
        print("INI file detected, start parsing...")
        config = ConfigParser.ConfigParser({"emin":"-3.0","emax":"3.0"})
        config.read(optFile)
        
        kpdata = {} # dict store the kpoints data
        for k, v in config.items("kpoints"):
            newk = tuple(map(lambda x: float(x), v.split()))
            newv = k.upper()
            if newk in kpdata: continue
            else: kpdata[newk] = newv
   
        print("Start reading vasprun.xml...") 
        xmldata = etree.parse("vasprun.xml")
        efermi = float(xmldata.xpath("//i[@name='efermi']/text()")[0])
        print("Fermi energy is %f eV" %efermi)
        tmp = xmldata.xpath("//varray[@name='kpointlist']/v")
        kplist = [] # preserve kpoints list
        for data in tmp:
            kplist.append(tuple(map(lambda x: float(x), data.text.split())))
        
        start = config.getint("segments","start")
        end = config.getint("segments","end")

        if end == -1: end = len(kplist) - 1 
        targetkp = kplist[start:end+1]

        print("K-POINTS read, in total %d" %len(targetkp))
        banddata = []
        for cnt in range(start,end+1):
            tmpdata = xmldata.xpath("//set[@comment='spin %d']/\
                                     set[@comment='kpoint %d']/r/text()" %(spin,cnt+1))
            # discard the data of occupation, and make the energy relative to Fermi energy
            banddata.append(map(lambda x: float(x.split()[0]) - efermi, tmpdata))
            #print tmpdata
        print("Band read, in total %d" %len(banddata[0]))

        emax = config.getfloat("erange","emax")
        emin = config.getfloat("erange","emin")

        print("Start plotting graph, from %.2f eV to %.2f eV" %(emin,emax))

        # with respect to the fermi energy
        maxline = 0 # set max to min
        minline = len(banddata[0]) - 1 # set min to max

        for data in banddata:
            fdata = filter(lambda x : x > emin and x < emax, data)
            maxn = data.index(max(fdata))
            minn = data.index(min(fdata))
            if maxn > maxline: maxline = maxn
            if minn < maxline: minline = minn
            # count the number of line needed
            #maxn = len(filter(lambda x: x>= 0, fdata))
            #minn = len(filter(lambda x: x< 0, fdata))
        print("Plotting from line %d to %d, in total %d lines" %(minline,maxline,maxline-minline+1))
        


        rec_basis = []
        for rec in xmldata.xpath("//structure[@name='finalpos']/crystal\
                                   /varray[@name='rec_basis']/v/text()"):
            #print rec
            rec_basis.append(map(lambda x: float(x), rec.split()))

        rec_basis = np.array(rec_basis).T
        #print "reciprocal vector", rec_basis
        
        currentX = 0.0
        kpstack = []
        tmpband = []

        cnt = 0
        xlabels = []

        #print banddata[0][minline:maxline+1]
        #print banddata[80][minline:maxline+1]
        for kp, ene in zip(targetkp,banddata):
            cnt += 1
            tmpband.append(ene[minline:maxline+1])

            if kp in kpdata:
                print("Find high symmetry points: %s" %kpdata[kp])
                #print cnt
                kpstack.append((kp,kpdata[kp]))
            
            Ytot = []
            if len(kpstack) == 2:
                startkp = kpstack[0]
                endkp = kpstack[1]
                print("Drawing lines from %s to %s" %(startkp[1],endkp[1]))

                vec = np.array(endkp[0]) - np.array(startkp[0])
                dist = np.linalg.norm(np.dot(rec_basis,vec))
                # print "distance",dist

                startX = currentX
                currentX += dist

                X = np.linspace(startX,currentX,len(tmpband),endpoint=True)

                #print len(tmpband)
                tmpband_T = zip(*tmpband)

                # not Ytot:
                #    for Y in tmpband_T:
                #        Ytot.append(copy.deepcopy(Y))
                #else:
                #    for ii, Y in enumerate(tmpband_T):
                #        Ytot[ii].extend(Y)
                #print tmpband_T[0]
                for Y in tmpband_T:
                    #newX, newY = smooth(X,Y)
                    #newY = sp.fft(Y)
                    #newY[2:] = 0.0
                    #Y = np.real(sp.ifft(newY))
                    #print(len(X),len(Y))
                    #func = interp1d(np.array(X),np.array(Y),kind='cubic')
                    #print X
                    #print Y[0]
                    #newX = np.linspace(min(X),max(X),9,endpoint=True)
                    #newY = func(newX)
                    plt.plot(X,Y,color='k')
                    
                if not xlabels:
                    xlabels = [(startX,startkp[1]),(currentX,endkp[1])]
                else:
                    xlabels.append((currentX,endkp[1]))

                #print xlabels
                kpstack = []
                tmpband = []

        #for Y in Ytot:
        #    print(len(Y))
        #    print len(X)
        #    func = interp1d(X,Y)
        #    newX = np.linspace(0,currentX,1000)
        #    newY = func(newX)
        #    plt.plot(newX,newY,color='k')
        for rec in xlabels:
            X = (rec[0],rec[0])
            Y = (emin,emax)
            plt.plot(X,Y,ls='--',color='k')
        
        plt.xticks(*zip(*xlabels))
        plt.xlim((0,currentX))
        plt.ylim((emin,emax))
        plt.show()        

    
def main_gui():
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Plot the band structure.",\
                                     prog = sys.argv[0])
    parser.add_argument('-m','--mode', default = 'cmd', type = str)
    parser.add_argument('-i','--ini', default = 'config.ini', type = str)

    res = parser.parse_args(sys.argv[1:])

    if res.mode == 'cmd':
        print("Using command mode...")
        if res.ini == parser.get_default('ini'):
            print("Notion: will use default init file: config.ini")
        main_cmd(optFile = res.ini)

    elif res.mode == 'gui':
        print("Using gui mode..., parameters obtained from gui windows")
        main_gui()
