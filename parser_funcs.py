import math
import os
import sys

#Units
gev=1;mev=1e-3*gev;kev=1e-6*gev;

#constants
hbar=float(1.054*1e-34/(1.6e-19)/(1e9));
speed_of_light=3e8;conversion=hbar**2*speed_of_light**2;
alpha_em=1.0/137.035999074;

#masses
mp=0.938272046*gev;melec=511*kev;mmuon=105.658*mev;
mpi=134.9767*mev;mpip=139.57018*mev;mkaon=0.493667*gev;
mj_psi=3.097*gev;

def dp(a1,a2):
    tot=0
    for i in range(len(a1)):
        tot+=a1[i]*a2[i]
    return tot
def sub(a1,a2):
    return [a1[i]-a2[i] for i in range(len(a1))]
def add(a1,a2):
    return [a1[i]+a2[i] for i in range(len(a1))]
def mul(c,a1):
    return [c*a1[i] for i in range(len(a1))]
def plane_cross(p,o,r,n):
    bot=dp(n,p)
    if bot == 0:
        return -1
    a = dp(n,sub(r,o))/bot
    R=sub(r,add(o,mul(a,p)))
    return dp(R,R)**0.5

def find_total(summary_file,run_num):
    with open(summary_file) as sf:
        dats=sf.read().splitlines()
        for i in range(len(dats)):
            line=dats[i].split()
            if len(line)==2 and line[1]==run_num:
                for j in range(len(dats)-i-1):
                    line2 = dats[i+j+1].split()
                    if len(line2) > 4 and line2[0]=="Total":
                        return float(line2[3])
        print("Could not find run " + run_num)
        return -1

def mom(arr):
    return math.sqrt(arr[3]**2+arr[1]**2+arr[2]**2)
    
def theta(arr):
    return math.acos(arr[3]/mom(arr))

def phi(arr):
    return math.atan(arr[0]/arr[1])

def cos_theta(arr):
    return arr[3]/mom(arr)

def theta_dif(arr1,arr2):
    return math.acos(dp(arr1[1:4],arr2[1:4])/mom(arr1)/mom(arr2))

#Expect each element of arr to be 4 elements. Expects px py pz E.
def invariant_mass(arr):
    tot = [0,0,0,0]
    for line in arr:
        for i in range(4):
            tot[i]+=line[i]
    return tot[0]**2-tot[3]**2-tot[1]**2-tot[2]**2

def cos_angle_between(arr1,arr2):
    return (arr1[3]*arr2[3]+arr1[1]*arr2[1]+arr1[2]*arr2[2])/mom(arr1)/mom(arr2)

def record_list(outfile,data):
    with open(outfile,'w') as w:
        for line in data:
            tmpstr=""
            for elem in line:
                tmpstr+=str(elem)+' '
            tmpstr=tmpstr[:-1]
            tmpstr+='\n'
            w.write(tmpstr)
def append_list(outfile,data):
    with open(outfile,'a') as w:
        for line in data:
            tmpstr=""
            for elem in line:
                tmpstr+=str(elem)+' '
            tmpstr=tmpstr[:-1]
            tmpstr+='\n'
            w.write(tmpstr)

#l_i = dist_i/speed_of_light/lifetime*mass
def dec_prob(l1,l2):
    return math.exp(-l1)-math.exp(-l2)
def dec_loc(l1,l2):
    try:
        return -math.log(math.exp(-l1)-random.random()*(math.exp(-l1)-math.exp(-l2)))
    except:
        return l1
def solve_dec_loc(l1,l2,pos):
    return (10**pos-math.exp(-l1))/(math.exp(-l2)-math.exp(-l1))
def dec_loc_set(l1,l2,pos):
    try:
        return -math.log(math.exp(-l1)-pos*(math.exp(-l1)-math.exp(-l2)))
    except:
        return l1
