#!/root/anaconda3/bin/python
from __future__ import print_function
import sys, os
import argparse
from itertools import islice
from subprocess import call
import numpy as np
import math
import re

# functions
#
def get_q1_val(step,idx,idxfix,dist): 
    # assume Q1 neg plus quadrant (left up) height[0,4]
    # vary width then height - follow slicing idea
    # initialize to zero
    hdown = 0.0 
    hup = 0.0 
    height = str(hdown)+","+str(hup)
    win = 0.0 
    wout = 0.0 
    width = str(win)+","+str(wout)
    # stop criterium
#   thresh  = float(0.0000050)  # height
#   thresh  = float(0.000050)  # height
    thresh  = float(0.000500)  # height
#   collect output to check if necessary space has been covered
    foutput = "jval_q1.txt"
    f1 = open(foutput,"w")

    #start Q1
    wlogic = True 
    i=1 
    valq1 = float(0.0)
    while (wlogic == True):
        w1 = float(win -i*step)
        wi = f"{w1:.3f}"
        w2 = float(wout - (i-1)*step)
        wo = f"{w2:.3f}"
        width = str(wi)+","+str(wo)
        l=1
        val = float(0.0)
        flag = True 
        while (flag == True):
            h1 = float(hdown + (l-1)*step)
            # store as string with precicion .3
            hd = f"{h1:.3f}"
            h2 = float(hup + l*step)
            hu = f"{h2:.3f}"
            height = str(hd)+","+str(hu)
#           print("height", height)
#           print("width", width)

            jval = get_jval(idx,idxfix,dist,height,width)
            sign3 = get_sign(jval)
#           print("jval", jval)
            # get center point of pixel
            xyz = find_center() 
#           print("xyz", xyz )
            # collect what we have to plot
            f1.write("%f %f %f %f \n" % (xyz[0],xyz[1],xyz[2],jval))

            ht = float(hdown + l*step)
            hdt = f"{ht:.3f}"
            htt = float(hup + (l+1)*step)
            hut = f"{htt:.3f}"
            htemp = str(hdt)+","+str(hut)
            jtmp = get_jval(idx,idxfix,dist,htemp,width)
            sign4 = get_sign(jtmp)
#           print("height", height)
#           print("sign 3", sign3)
#           print("htemp", htemp)
#           print("sign 4", sign4)
#           print("l", l)
#           assume step 0.2 then 1 Angstroem buffer
            if (l > 13) and (float(sign4) != float(sign3)):
                print("inner loop if l block")
                flag = False
                break 
            elif ((abs(float(jval))) <= thresh):
                print("inner loop if jval block")
                flag = False
                break
            elif (l > 28):
                print("q1 l larger 28 exit loop")
                flag = False
                break
            # sum what we have
            val = val + float(jval) 
#           print("val", val)
            l = l + 1 
    
        valq1 = valq1 + val 
#       print("val Q1 ", valq1) 
        #
        jval1 = get_jval(idx,idxfix,dist,"0.0,0.1",width)
        sign1 = get_sign(jval1)
#       print("after jval1" )
        t1 = float(win-i*step)  
        ti = f"{t1:.3f}"
        t2 = float(wout-(i+1)*step)
        to = f"{t2:.3f}"
        temp = str(ti)+","+str(to)
#       temp = str(win-(i-1)*step)+","+str(wout-i*step)
        jval2 = get_jval(idx,idxfix,dist,"0.0,0.1",temp)
        sign2 = get_sign(jval2)
#       print("after jval2" )
#       print("sign1", sign1) 
#       print("sign2", sign2) 
#       print("i", i) 
#       print("l", l) 

        # stop criterium for width
        if (i > 10) and (float(sign2) != float(sign1)):
            print("I am in outer loop if sign block")
            wlogic = False
            break 
        elif (abs(val) <= thresh):
            print("I am in outer loop if abs val block")
            wlogic = False
            break 
        elif (i > 20):
            print("q1 i lager 20 exit loop")
            wlogic = False
            break
        # counter
        i = i + 1

    print("valq1 final =", valq1) 
    f1.close()
    return valq1
#
def get_q2_val(step,idx,idxfix,dist): 
    # assume Q2 neg neg quadrant (left down) height[-4,0]
    # vary width then height - follow slicing idea
    # initialize to zero
    hdown = 0.0 
    hup = 0.0 
    height = str(hdown)+","+str(hup)
    win = 0.0 
    wout = 0.0 
    width = str(win)+","+str(wout)
    # stop criterium
#   thresh  = float(0.0000050)  # height
#   thresh  = float(0.000050)  # height
    thresh  = float(0.000500)  # height
#   collect output to check if necessary space has been covered
    foutput = "jval_q2.txt"
    f1 = open(foutput,"w")

    #start Q2
    wlogic = True 
    i=1 
    valq2 = float(0.0)
    while (wlogic == True):
        w1 = float(win -i*step)
        wi = f"{w1:.3f}"
        w2 = float(wout - (i-1)*step)
        wo = f"{w2:.3f}"
        width = str(wi)+","+str(wo)
        l=1
        val = float(0.0)
        flag = True 
        while (flag == True):
         #   height = str(hdown-(l-1)*step)+","+str(hup-l*step)
            h1 = float(hdown - (l-1)*step)
            # store as strin with precicion .3
            hd = f"{h1:.3f}"
            h2 = float(hup - l*step)
            hu = f"{h2:.3f}"
            height = str(hu)+","+str(hd)
#           height = str(hd)+","+str(hu)
#           print("height", height)
#           print("width", width)

            jval = get_jval(idx,idxfix,dist,height,width)
            sign3 = get_sign(jval)
#           print("jval", jval)
            # get center point of pixel
            xyz = find_center() 
#           print("xyz", xyz )
            # collect what we have to plot
            f1.write("%f %f %f %f \n" % (xyz[0],xyz[1],xyz[2],jval))

            ht = float(hdown - l*step)
            hdt = f"{ht:.3f}"
            htt = float(hup - (l+1)*step)
            hut = f"{htt:.3f}"
            htemp = str(hdt)+","+str(hut)
            jtmp = get_jval(idx,idxfix,dist,htemp,width)
            sign4 = get_sign(jtmp)
#           print("height", height)
#           print("width", width)
#           print("sign 3", sign3)
#           print("htemp", htemp)
#           print("sign 4", sign4)
#           print("l", l)
#           print("i", i)
#           assume step 0.2 then 1 Angstroem buffer
            if (l > 13) and (float(sign4) != float(sign3)):
                print("inner loop if l block")
                flag = False
                break 
            elif ((abs(float(jval))) <= thresh):
                print("inner loop if jval block")
                flag = False
                break
            elif (l > 28):
                print("q2 l larger 28 exit loop")
                flag = False
                break 

            # sum what we have
            val = val + float(jval) 
#           print("val", val)
            l = l + 1 
    
        valq2 = valq2 + val 
#       print("val Q2 ", valq2) 
        #
        jval1 = get_jval(idx,idxfix,dist,"-0.1,0.0",width)
        sign1 = get_sign(jval1)

        t1 = float(win- i*step)  
        ti = f"{t1:.3f}"
        t2 = float(wout-(i+1)*step)
        to = f"{t2:.3f}"
        temp = str(ti)+","+str(to)
#       temp = str(win-(i-1)*step)+","+str(wout-i*step)
        jval2 = get_jval(idx,idxfix,dist,"-0.1,0.0",temp)
        sign2 = get_sign(jval2)
#       print("after jval2" )
#       print("val", val )
#       print("sign1", sign1) 
#       print("sign2", sign2) 
#       print("i", i)
#       print("wlogic outer loop", wlogic)
        # stop criterium for width is sign change
        if (i > 10) and (float(sign2) != float(sign1)):
            print("outer loop if sign block")
            wlogic = False
            break 
        elif ((abs(float(val))) <= thresh):
            print("outer loop if thresh block")
            wlogic = False
            break 
        elif (i > 20):
            print("q2 i larger 20 exit loop")
            wlogic = False
            break 
#       print("wlogic end outer loop", wlogic)
#       counter
        i = i + 1

    print("valq2 final =", valq2) 
    f1.close()
    return valq2

def get_q3_val(step,idx,idxfix,dist): 
    # assume Q3 plus plus quadrant (right up) height[0,4], width[0,4]
    # vary width then height - follow slicing idea
    # initialize to zero
    hdown = 0.0 
    hup = 0.0 
    height = str(hdown)+","+str(hup)
    win = 0.0 
    wout = 0.0 
    width = str(win)+","+str(wout)
    # stop criterium
#   thresh  = float(0.0000050)  # height
#   thresh  = float(0.000050)  # height
    thresh  = float(0.000500)  # height
#   collect output to check if necessary space has been covered
    foutput = "jval_q3.txt"
    f1 = open(foutput,"w")

    #start Q3
    wlogic = True 
    i=1 
    valq3 = float(0.0)
    while (wlogic == True):
        w1 = float(win +i*step)
        wi = f"{w1:.3f}"
        w2 = float(wout + (i-1)*step)
        wo = f"{w2:.3f}"
        width = str(wo)+","+str(wi)
        l=1
        val = float(0.0)
        flag = True 
        while (flag == True):
            h1 = float(hdown + (l-1)*step)
            # store as strin with precicion .3
            hd = f"{h1:.3f}"
            h2 = float(hup + l*step)
            hu = f"{h2:.3f}"
            height = str(hd)+","+str(hu)
#           print("height", height)
#           print("width", width)

            jval = get_jval(idx,idxfix,dist,height,width)
            sign3 = get_sign(jval)
#           print("jval", jval)
            # get center point of pixel
            xyz = find_center() 
#           print("xyz", xyz )
            # collect what we have to plot
            f1.write("%f %f %f %f \n" % (xyz[0],xyz[1],xyz[2],jval))

            ht = float(hdown + l*step)
            hdt = f"{ht:.3f}"
            htt = float(hup + (l+1)*step)
            hut = f"{htt:.3f}"
            htemp = str(hdt)+","+str(hut)
            jtmp = get_jval(idx,idxfix,dist,htemp,width)
            sign4 = get_sign(jtmp)
#           print("height", height)
#           print("sign 3", sign3)
#           print("htemp", htemp)
#           print("sign 4", sign4)
#           print("l", l)
#           assume step 0.2 then 1 Angstroem buffer
            if (l > 13) and (float(sign4) != float(sign3)):
#               print("inner loop if l block")
                flag = False
                break 
            elif ((abs(float(jval))) <= thresh):
#               print("inner loop if jval block")
                flag = False
                break
            elif (l > 28):
                print("q3 l larger 28 exit loop")
                flag = False
                break

            # sum what we have
            val = val + float(jval) 
#           print("val", val)
            l = l + 1 
    
        valq3 = valq3 + val 
#       print("val Q3 ", valq3) 
        #
#       print("width ", width)
        jval1 = get_jval(idx,idxfix,dist,"0.0,0.1",width)
        sign1 = get_sign(jval1)
#       print("sign1", sign1)
        t1 = float(wout + i*step)  
        ti = f"{t1:.3f}"
        t2 = float(win + (i+1)*step)
        to = f"{t2:.3f}"
        temp = str(ti)+","+str(to)
#       temp = str(wout+ (i-1)*step)+","+str(win + i*step)
#       print("temp ", temp)
        jval2 = get_jval(idx,idxfix,dist,"0.0,0.1",temp)
        sign2 = get_sign(jval2)
#       print("sign2", sign2)
        if (i > 10) and (float(sign2) != float(sign1)):
            wlogic = False
            break 
        elif (abs(val) <= thresh):
            wlogic = False
            break 
        elif (i > 20):
            print("q3 i larger 20 exit loop")
            wlogic = False
            break
        # counter
        i = i + 1

    print("valq3 final =", valq3) 
    f1.close()
    return valq3
#
def get_q4_val(step,idx,idxfix,dist): 
    # assume Q4 plus minus quadrant (right down) height[-4,0], width[0,4]
    # vary width then height - follow slicing idea
    # initialize to zero
    hdown = 0.0 
    hup = 0.0 
    height = str(hdown)+","+str(hup)
    win = 0.0 
    wout = 0.0 
    width = str(win)+","+str(wout)
    # stop criterium
#   thresh  = float(0.0000050)  # height
#   thresh  = float(0.000050)  # height
    thresh  = float(0.000500)  # height
#   collect output to check if necessary space has been covered
    foutput = "jval_q4.txt"
    f1 = open(foutput,"w")

    #start Q4
    wlogic = True 
    i=1 
    valq4 = float(0.0)
    while (wlogic == True):
        w1 = float(win +i*step)
        wi = f"{w1:.3f}"
        w2 = float(wout + (i-1)*step)
        wo = f"{w2:.3f}"
        width = str(wo)+","+str(wi)
        l=1
        val = float(0.0)
        flag = True 
        while (flag == True):
            h1 = float(hdown - (l-1)*step)
            # store as strin with precicion .3
            hd = f"{h1:.3f}"
            h2 = float(hup - l*step)
            hu = f"{h2:.3f}"
            height = str(hu)+","+str(hd)
#           print("height", height)
#           print("width", width)

            jval = get_jval(idx,idxfix,dist,height,width)
            sign3 = get_sign(jval)
#           print("jval", jval)
            # get center point of pixel
            xyz = find_center() 
#           print("xyz", xyz )
            # collect what we have to plot
            f1.write("%f %f %f %f \n" % (xyz[0],xyz[1],xyz[2],jval))

            ht = float(hdown - l*step)
            hdt = f"{ht:.3f}"
            htt = float(hup - (l+1)*step)
            hut = f"{htt:.3f}"
            htemp = str(hdt)+","+str(hut)
            jtmp = get_jval(idx,idxfix,dist,htemp,width)
            sign4 = get_sign(jtmp)
#           print("height", height)
#           print("sign 3", sign3)
#           print("htemp", htemp)
#           print("sign 4", sign4)
#           print("l", l)
#           assume step 0.2 then 1 Angstroem buffer
            if (l > 13) and (float(sign4) != float(sign3)):
#               print("inner loop if l block")
                flag = False
                break 
            elif ((abs(float(jval))) <= thresh):
#               print("inner loop if jval block")
                flag = False
                break
            elif (l > 28):
                print("q4 l larger 28 exit loop")
                flag = False
                break
            # sum what we have
            val = val + float(jval) 
#           print("val", val)
            l = l + 1 
    
        valq4 = valq4 + val 
#       print("val Q4 ", valq4) 
        #
        jval1 = get_jval(idx,idxfix,dist,"-0.1,0.0",width)
        sign1 = get_sign(jval1)
        t1 = float(wout + i*step)  
        ti = f"{t1:.3f}"
        t2 = float(win + (i+1)*step)
        to = f"{t2:.3f}"
        temp = str(ti)+","+str(to)
#       temp = str(wout+ (i-1)*step)+","+str(win + i*step)
        jval2 = get_jval(idx,idxfix,dist,"-0.1,0.0",temp)
        sign2 = get_sign(jval2)
        # stop criterium for width
        if (i > 10) and (float(sign2) != float(sign1)):
            wlogic = False
            break 
        elif (abs(val) <= thresh):
            wlogic = False
            break 
        elif (i > 20):
            print("q4 i larger 20 exit loop")
            wlogic = False
            break
        # counter
        i = i + 1

    print("valq4 final =", valq4) 
    f1.close()
    return valq4
#


def get_sign(jval):
    val = float(jval)
    if (val < 0.0):
        sign = -1.0  
    if (val > 0.0):
        sign = +1.0
    if (val == 0.0):
        sign = "none"
    return sign 

def get_jval(idx,idxfix,dist,height,width):
    # ouput j as float
    # set up gimic input
    finput = "gimic.inp"
    # height = "-3.0,3.0"
    # width = "-2.0,2.0" 
    write_gimic_input(finput,idx,idxfix,dist,height,width)
    # run gimic
    os.system("gimic > gimic.out ")
#   os.system("~/source/gimic/build/bin/gimic > gimic.out")
    # get output
    jval = find_j()
#   print(jval)
    return jval

def find_j():
    pattern = "Induced current (nA/T)  : "
    left = "Left handed" 
    filename = "gimic.out" 
    fac=1.0
    with open(filename) as f:
        for line in f:
            line = line.rstrip() 
            if left in line:
                fac = -1.0 
                print("left handed") 
            if pattern in line:
#               print(line) 
                l = line.strip().split()
                j = float(l[4])*float(fac) 
                return j 

def find_center():
    fin = "grid.xyz"
    pattern = "X"
    xcoord = []
    ycoord = []
    zcoord = []
    xyz = []
    with open(fin) as f:
        for line in f:
            line = line.rstrip()
            if pattern in line:
                l = line.strip().split()
                xcoord.append(float(l[1]))
                ycoord.append(float(l[2]))
                zcoord.append(float(l[3]))
    xsum = float(0.0)
    ysum = float(0.0)
    zsum = float(0.0)
    for i in range(3):
        xsum = xsum + xcoord[i]
        ysum = ysum + ycoord[i]
        zsum = zsum + zcoord[i]
    
    xyz.append(xsum/4.0)
    xyz.append(ysum/4.0)
    xyz.append(zsum/4.0)
    # print("xyz", xyz) 
    return xyz


def write_gimic_input(finput,idx,idxfix,dist,height,width):
    f1 = open(finput,'w')
    f1.write("# GIMIC INPUT \n")
    f1.write("calc=integral \n")
    f1.write('title="" \n')
    f1.write('basis="MOL" \n')
    f1.write('xdens="XDENS" \n')
    f1.write("debug=1 \n")
    f1.write("openshell=false \n")
    f1.write("magnet=[0.0,0.0,1.0] \n")
#   f1.write("magnet_axis=X \n")
    f1.write("Grid(bond) { \n")
    f1.write("    type=gauss \n")
    # magnet_axis=X 

    idx = map(str,idx)
    line = ",".join(idx)
    string = '    bond=[',line,']' 
    string = map(str,string)
    newline = "".join(string) 
    f1.write(newline + "\n")

    string = '    fixpoint=',str(idxfix) 
    newline = "".join(string) 
    f1.write(newline +"\n")
    f1.write("    distance=" + str(dist) + "\n")
    f1.write("    gauss_order=9 \n")
    f1.write("    spacing=[0.01, 0.01, 0.01] \n")

    height = map(str,height)
    line = "".join(height)
    string = '    height=[',line,']'
    newline = "".join(string) 
    f1.write(newline + "\n")

    width = map(str,width)
    line = "".join(width)
    string = '    width=[',line,']'
    newline = "".join(string) 
    f1.write(newline + "\n")
    f1.write("} \n")
    f1.write(" \n")
    f1.write("Advanced { \n") 
    f1.write("    lip_order=5 \n") 
    f1.write("    spherical=off \n") 
    f1.write("    diamag=on \n") 
    f1.write("    paramag=on \n") 
    f1.write("    GIAO=on  \n") 
    f1.write("    screening=on \n") 
    f1.write("    screening_thrs=1.d-8  \n") 
    f1.write("} \n") 
    f1.close()

def get_dist(idx):
    ang2bohr = 1.8897161646320724
    bohr2ang = 0.52918
    fin = "coord.xyz"
    #
    natom = 0
    element = []
    coord_x = []
    coord_y = []
    coord_z = []
    with open(fin) as f:
        # skip first two lines
        for _ in range(2):
            next(f)
        for i in islice(f, 0, None):
            natom = natom + 1 
            l = i.strip().split()
            element.append(l[0])
            coord_x.append(float(l[1]))
            coord_y.append(float(l[2]))
            coord_z.append(float(l[3]))
    # fin.close()
    print("Number of atoms:", natom)
    print("Bond defined by atoms: ", idx[0], idx[1])

    # A
    # print("coord atom",idx[0],": ", coord_x[idx[0]-1], coord_y[idx[0]-1], coord_z[idx[0]-1])
    # B
    #  print("coord atom",idx[1],": ", coord_x[idx[1]-1], coord_y[idx[1]-1], coord_z[idx[1]-1])
    # get bond vector  B-A = vec A->B
    dx = coord_x[idx[1]-1] - coord_x[idx[0]-1]
    dy = coord_y[idx[1]-1] - coord_y[idx[0]-1]
    dz = coord_z[idx[1]-1] - coord_z[idx[0]-1]
    # get bond distance
    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    print("bond distance in A", dist)
    print("bond distance in bohr", dist*ang2bohr)
    # save 0.5 dist as gimic input in bohr
    halfdist_bohr = 0.5*dist*ang2bohr 
    return halfdist_bohr

