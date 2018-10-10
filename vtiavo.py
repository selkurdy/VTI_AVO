"""Model AVO incorporating VTI Thomsen parameters expsilon and delta."""

# -*- coding: utf-8 -*-
import sys, os.path
import os
import datetime
import argparse
import numpy as np
import numpy.ma as ma
import math as m
import matplotlib.pyplot as plt

"""
Created on Fri Mar 06 21:37:11 2015

@author: Sami Elkurdy

AVO modeling/Inversion with VTI
Mar 08 2015 added graphing for both iso and vti
    epsilons is the SV wave anisotropy: deviation of phase velocity surface from ellipticity.
    For elliptically anisotropic media espsilons = 0

Mar 11 2015 added phase velocity computation vp and vs
+ Apr 7 2015
    added option to enter delta (from hampson Russell) or use epsilonS to model velocity with angle
    modified plotting to have 2 subplots
+ to enter epsilon s
>python vtiavo.py --model vva_epes --vp 2000 3000 --vs 1200 1500  --delta 0.2 0.1  --epsilonp 0.2 0.22

+ to enter delta
>python vtiavo.py --model vva_delta --vp 2000 3000 --vs 1200 1500  --delta 0.2 0.1  --epsilonp 0.2 0.22
>python vtiavo.py --model vti --vp 2000 3000 --vs 1200 1500  --delta 0.2 0.1  --epsilonp 0.2 0.22
>python vtiavo.py --model iso --vp 2000 3000 --vs 1200 1500  --delta 0.2 0.1  --epsilonp 0.2 0.22

"""
def d2r(ang):
    """convert degrees to radians."""
    angr = (ang * np.pi / 180.0)
    return angr


def r2d(ang):
    """convert radians to degrees."""
    return (ang * 180.0 / np.pi)

# single angle computation for later use of inversion
def Ripp1(vp1,vp2,vs1,vs2,rb1,rb2,ang):
    """Modified version of reflectivity with angle."""
    dvp = vp2 - vp1
    avp = (vp2 + vp1) / 2.0
    dvs = vs2 - vs1
    avs = (vs2 + vs1) / 2.0
    drb = rb2 - rb1
    arb = (rb2 + rb1) / 2.0
    angr = d2r(ang)
    vsvpsqr = (avs / avp) ** 2
    ri1 = 0.5 * (1.0 - 4.0 * vsvpsqr * m.sin(angr)**2) * drb / arb
    ri2 = (1.0 / m.cos(angr)**2) * dvp / avp
    ri3 = 4.0 * vsvpsqr * m.sin(angr)**2 * dvs / avs
    return (ri1 + ri2 - ri3)


def Ripp(vp1,vp2,vs1,vs2,rb1,rb2,ang):
    """Compute reflectivity with angle."""
    dvp = vp2 - vp1
    avp = (vp2 + vp1) / 2.0
    dvs = vs2 - vs1
    avs = (vs2 + vs1) / 2.0
    drb = rb2 - rb1
    arb = (rb2 + rb1) / 2.0
    angr = d2r(ang)
    vsvpsqr = (avs / avp) ** 2
    ri1 = 0.5 * (1.0 - 4.0 * vsvpsqr * np.sin(angr)**2) * drb/arb
    ri2 = (1.0 / np.cos(angr)**2) * dvp /avp
    ri3 = 4.0 * vsvpsqr * np.sin(angr)**2 * dvs/avs
    return (ri1 + ri2 - ri3)


def Rapp(vp1,vp2,vs1,vs2,rb1,rb2,dd1,dd2,ang):
    """Modified version of reflectivity with angle."""
    dvp = vp2 - vp1
    avp = (vp2 + vp1)/ 2.0
    dvs = vs2 - vs1
    avs = (vs2 + vs1) / 2.0
    drb = rb2 - rb1
    arb = (rb2 + rb1) / 2.0
    angr = d2r(ang)
    ddd = dd2 - dd1
    vsvpsqr = (avs / avp) ** 2
    ri1 = 0.5 * (1.0 - 4.0 * vsvpsqr * np.sin(angr)**2) * drb / arb
    ri2 = (1.0 / np.cos(angr)**2) * dvp / avp
    ri3 = 4.0 * vsvpsqr * np.sin(angr)**2 * dvs / avs
    ra = ddd / 2.0 * np.sin(angr)**2
    ri = ri1 + ri2 - ri3
    return ri,ri+ra

def Rapp1(vp1,vp2,vs1,vs2,rb1,rb2,dd1,dd2,ang):
    """single angle computation for later use of inversio."""
    dvp = vp2 - vp1
    avp = (vp2 + vp1) / 2.0
    dvs = vs2 - vs1
    avs = (vs2 + vs1) / 2.0
    drb = rb2 - rb1
    arb = (rb2 + rb1) / 2.0
    angr = d2r(ang)
    ddd = dd2 - dd1
    vsvpsqr = (avs / avp) ** 2
    ri1 = 0.5 * (1.0 - 4.0 * vsvpsqr * m.sin(angr)**2) * drb / arb
    ri2 = (1.0 / m.cos(angr)**2) * dvp / avp
    ri3 = 4.0 * vsvpsqr * m.sin(angr)**2 * dvs / avs
    ra = ddd / 2.0 * m.sin(angr)**2
    ri = ri1 + ri2 - ri3
    return ri + ra
    # returns only one value of anisotropic amplitude at a specific angle

def vapp(vp0,vs0,ep,es,ang):
    """equation taken from "Seismic Reflection Processing" Springer."""
    angr = d2r(ang)
    vsvpsqr = (vs0/vp0) ** 2
    vapp0= ep * np.sin(angr)**4
    vapp1= (ep - 4.0 * vsvpsqr) * np.sin(angr)**2 * np.cos(angr)**2
    vapp2 = vp0 * ( 1.0 + vapp1 + vapp0)
    return vapp2

def vapp_delta(vp0,ep,delta,ang):
    """equation taken from HampsonRussell trainaing manual."""
    angr = d2r(ang)
    vapp0 = ep * np.sin(angr)**4
    vapp1 = delta * np.sin(angr)**2 * np.cos(angr)**2
    vapp2 = vp0 * (1.0 + vapp1 + vapp0)
    return vapp2

def vasv(vs0,es,ang):
    """equation taken from "Seismic Reflection Processing" Springer."""
    angr = d2r(ang)
    return vs0 *(1.0 + 4.0 * es * np.sin(angr)**2 *np.cos(angr)**2)

def vasv_delta(vp0,vsv0,ep,delta,ang):
    """equation taken from HampsonRussell trainaing manual."""
    angr = d2r(ang)
    vpvssqr = (vp0/vsv0) ** 2
    vsvtheta = vsv0 * (1.0 + vpvssqr * (ep - delta) * np.sin(angr)**2 * np.cos(angr)**2)
    return vsvtheta

def delta(vp0,vs0,ep,es):
    """equation taken from HampsonRussell trainaing manual."""
    return (ep - 4.0 * es * (vs0 / vp0)**2)

def list_ri(a,r):
    """list angle with isotropic amplitude."""
    for i in range(a.size):
        print("%10.4f  %10.4f" % (a[i],r[i]))

def list_ra(a,ri,ra):
    """list angle vs isotropic and vti amplitudes."""
    for i in range(a.size):
        print("%10.4f  %10.4f  %10.4f" % (a[i],ri[i],ra[i]))

def getcommandline():
    """Get command line options."""
    vs1 = 4700.0
    vs2 = 5200.0
    vp1 = 9500.0
    vp2 = 10000.0
    rb1 = 2.70
    rb2 = 2.65
    amin = 0.0
    amax = 60.0
    epsilonp1 = 0.1
    epsilonp2 = 0.2
    epsilons1 = 0.05
    epsilons2 = 0.08
    delta1 = 0.1
    delta2 = 0.2

    parser = argparse.ArgumentParser(description='AVO w/ VTI modeling and inversion Mar 6, 2015 ')
    parser.add_argument('--vp',type=float,nargs=2,default=[vp1,vp2],help='Vp top bottom. dfv= %6.0f %6.0f' % (vp1,vp2))
    parser.add_argument('--vs',type=float,nargs=2,default=[vs1,vs2],help='Vs top bottom. dfv=%5.0f %5.0f' % (vs1,vs2))
    parser.add_argument('--rb',type=float,nargs=2,default=[rb1,rb2],help='Rhob top bottom. dfv=%5.2f %5.2f' % (rb1,rb2))
    parser.add_argument('--angles',type=float,nargs=2,default=[amin,amax],
        help='Minimum Maximum angles.dfv= %3.0f  %3.0f deg' % (amin,amax))
    parser.add_argument('--epsilonp',type=float, nargs=2,default=[epsilonp1,epsilonp2],
        help='Epsilon P of upper and lower layer. dfv=%5.2f  %5.2f' % (epsilonp1,epsilonp2))
    parser.add_argument('--epsilons',type=float, nargs=2,default=[epsilons1,epsilons2],
        help='Epsilon S of upper and lower layer. dfv=%5.2f  %5.2f' % (epsilons1,epsilons2))
    parser.add_argument('--delta',type=float, nargs=2,default=[delta1,delta2],
        help='Delta Delta of upper and lower layer. dfv=%5.2f  %5.2f' % (delta1,delta2))
    parser.add_argument('--model',choices=['iso','vti','vva_delta','vva_epsilon'],default='iso',
        help='modeling options:isotropic,anisotropic VTI, Phase velocities with angle epsilon, delta,\
        Phase velocities with angle epsilonP epsilonS -> delta is calculated .dfv= iso')
    parser.add_argument('--hideplot',action='store_true',default=False,
        help='Do not display plots, only save to pdf')

    results = parser.parse_args()
    return results

def main():
    """Main program."""
    cmdl = getcommandline()
    pdfiso0 = f'{cmdl.model}' + '_isoavo.pdf'
    pdfiso1 = f'{cmdl.model}' + '_igavo.pdf'
    pdfaniso0 = f'{cmdl.model}' + '_anavo.pdf'
    pdfaniso1 = f'{cmdl.model}' + '_anigavo.pdf'
    pdfvve = f'{cmdl.model}' + '_vveps_avo.pdf'
    pdfvvd = f'{cmdl.model}' + '_vvd_avo.pdf'

    anglearray = np.linspace(cmdl.angles[0],cmdl.angles[1],30,endpoint=True)
    anglearraysqr = np.power(np.sin(np.deg2rad(anglearray)),2)
    # print(anglearraysqr.size,anglearray.size)

    if cmdl.model == 'iso':
        ripp = Ripp(cmdl.vp[0],cmdl.vp[1],cmdl.vs[0],cmdl.vs[1],cmdl.rb[0],cmdl.rb[1],anglearray)
        list_ri(anglearray,ripp)
        plt.plot(anglearray,ripp,lw=3,label='ISO')
        plt.title('ISOTROPIC AVO')
        plt.annotate('Vp1:%6.0f Vs1:%6.0f Rb1:%3.2f' % (cmdl.vp[0],cmdl.vs[0],cmdl.rb[0]),
            xy=(anglearray[5],ripp[20]),xytext=(0.5,0.81),textcoords='figure fraction')
        plt.annotate('Vp2:%6.0f Vs2:%6.0f Rb2:%3.2f' % (cmdl.vp[1],cmdl.vs[1],cmdl.rb[1]),
            xy=(anglearray[5],ripp[20]),xytext=(0.5,0.75),textcoords='figure fraction')
        plt.ylabel('Reflection Amplitude')
        plt.xlabel('Incidence Angle in degrees')
        plt.legend(loc='lower left')
        # plt.legend(loc='lower center')
        plt.savefig(pdfiso0)
        if not cmdl.hideplot:
            plt.show()
        plt.plot(anglearraysqr,ripp,lw=3)
        plt.title('ISOTROPIC AVO SIN SQUARED INCIDENCE')
        plt.annotate('Vp1:%6.0f Vs1:%6.0f Rb1:%3.2f' % (cmdl.vp[0],cmdl.vs[0],cmdl.rb[0]),
            xy=(anglearraysqr[5],ripp[20]),xytext=(0.5,0.81),textcoords='figure fraction')
        plt.annotate('Vp2:%6.0f Vs2:%6.0f Rb2:%3.2f' % (cmdl.vp[1],cmdl.vs[1],cmdl.rb[1]),
            xy=(anglearraysqr[5],ripp[20]),xytext=(0.5,0.75),textcoords='figure fraction')
        plt.ylabel('Reflection Amplitude')
        plt.xlabel('Incidence Angle in radians squared')
        xi = np.linspace(anglearraysqr.min(),anglearraysqr.max())
        acf = np.polyfit(anglearraysqr,ripp,1)
        ampii = np.polyval(acf,xi)
        plt.plot(xi,ampii,c='r',lw=2,label='I_G')
        plt.annotate('Intercept:%6.3f Gradient:%6.3f ' % (acf[1],acf[0]),
            xy=(anglearraysqr[5],ripp[20]),xytext=(0.5,0.68),textcoords='figure fraction')
        plt.legend(loc='lower left')
        # plt.legend(loc='lower center')
        plt.savefig(pdfiso1)
        if not cmdl.hideplot:
            plt.show()

    if cmdl.model == 'vva_delta':
        vp1ang = vapp_delta(cmdl.vp[0],cmdl.epsilonp[0],cmdl.delta[0],anglearray)
        vp2ang = vapp_delta(cmdl.vp[1],cmdl.epsilonp[1],cmdl.delta[1],anglearray)
        vs1ang = vasv_delta(cmdl.vp[0],cmdl.vs[0],cmdl.epsilons[0],cmdl.delta[0],anglearray)
        vs2ang = vasv_delta(cmdl.vp[1],cmdl.vs[1],cmdl.epsilons[1],cmdl.delta[1],anglearray)

        plt.figure(1,figsize=(7,7))
        plt.subplot(2,1,1)
        plt.plot(anglearray,vp1ang,'b',label='Vp1',lw=3)
        plt.plot(anglearray,vp2ang,'g',label='Vp2',lw=3)
        plt.xlabel('Angle of Incidence')
        plt.title('P Phase Velocity vs Angle')
        plt.legend(loc='center right')
        plt.annotate('Vp1:%6.0f Vs1:%6.0f epsp1:%3.2f epss1:%3.2f d1:%3.2f' %
            (cmdl.vp[0],cmdl.vs[0],cmdl.epsilonp[0],cmdl.epsilons[0],cmdl.delta[0]),
            xy=(anglearray[5],vp1ang[20]),xytext=(0.3,0.70),textcoords='figure fraction')
        plt.annotate('Vp2:%6.0f Vs2:%6.0f epsp2:%3.2f epss2:%3.2f d2:%3.2f' %
            (cmdl.vp[1],cmdl.vs[1],cmdl.epsilonp[1],cmdl.epsilons[1],cmdl.delta[1]),
            xy=(anglearray[5],vp1ang[20]),xytext=(0.3,0.67),textcoords='figure fraction')
        plt.ylabel('Phase P Velocity')
        plt.subplot(2,1,2)
        plt.plot(anglearray,vs1ang,'r',label='Vs1',lw=3)
        plt.plot(anglearray,vs2ang,'m',label='Vs2',lw=3)
        plt.xlabel('Angle of Incidence')
        plt.ylabel('Phase SV Velocity')
        plt.legend(loc='center right')
        plt.title('SV Phase Velocity vs Angle')
        plt.annotate('Vp1:%6.0f Vs1:%6.0f epsp1:%3.2f epss1:%3.2f d1: %3.2f' %
            (cmdl.vp[0],cmdl.vs[0],cmdl.epsilonp[0],cmdl.epsilons[0],cmdl.delta[0]),
            xy=(anglearray[5],vs1ang[20]),xytext=(0.3,0.30),textcoords='figure fraction')
        plt.annotate('Vp2:%6.0f Vs2:%6.0f epsp2:%3.2f epss2:%3.2f d2:%3.2f' %
            (cmdl.vp[1],cmdl.vs[1],cmdl.epsilonp[1],cmdl.epsilons[1],cmdl.delta[1]),
            xy=(anglearray[5],vs1ang[20]),xytext=(0.3,0.27),textcoords='figure fraction')
        plt.tight_layout()
        plt.savefig(pdfvvd)
        if not cmdl.hideplot:
            plt.show()

    if cmdl.model == 'vva_epsilon':
        # equation taken from Seismic Data Processing, Springer Upaday
        vp1ang = vapp(cmdl.vp[0],cmdl.vs[0],cmdl.epsilonp[0],cmdl.epsilons[0],anglearray)
        vp2ang = vapp(cmdl.vp[1],cmdl.vs[1],cmdl.epsilonp[1],cmdl.epsilons[1],anglearray)
        vs1ang = vasv(cmdl.vs[0],cmdl.epsilons[0],anglearray)
        vs2ang = vasv(cmdl.vs[1],cmdl.epsilons[1],anglearray)

        plt.figure(1,figsize=(7,7))
        # plt.figure(1)
        plt.subplot(2,1,1)
        # plt.tight_layout()
        plt.plot(anglearray,vp1ang,'b',label='Vp1',lw=3)
        plt.plot(anglearray,vp2ang,'g',label='Vp2',lw=3)
        plt.xlabel('Angle of Incidence')
        plt.title('P Phase Velocity vs Angle')
        plt.legend(loc='center right')
        plt.annotate('Vp1:%6.0f Vs1:%6.0f epsp1: %3.2f, epss1: %3.2f' % (cmdl.vp[0],cmdl.vs[0],cmdl.epsilonp[0],cmdl.epsilons[0]),
            xy=(anglearray[5],vp1ang[20]),xytext=(0.3,0.70),textcoords='figure fraction')
        plt.annotate('Vp2:%6.0f Vs2:%6.0f  epsp2:%3.2f  epss2:%3.2f' % (cmdl.vp[1],cmdl.vs[1],cmdl.epsilonp[1],cmdl.epsilons[1]),
            xy=(anglearray[5],vp1ang[20]),xytext=(0.3,0.67),textcoords='figure fraction')
        plt.ylabel('Phase P Velocity')
        plt.subplot(2,1,2)
        plt.plot(anglearray,vs1ang,'r',label='Vs1',lw=3)
        plt.plot(anglearray,vs2ang,'m',label='Vs2',lw=3)
        plt.xlabel('Angle of Incidence')
        plt.ylabel('Phase SV Velocity')
        plt.legend(loc='center right')
        plt.title('SV Phase Velocity vs Angle')
        plt.annotate('Vp1:%6.0f Vs1:%6.0f epsp1: %3.2f epss1:%3.2f' %
            (cmdl.vp[0],cmdl.vs[0],cmdl.epsilonp[0],cmdl.epsilons[0]),
            xy=(anglearray[5],vs1ang[20]),xytext=(0.3,0.30),textcoords='figure fraction')
        plt.annotate('Vp2:%6.0f Vs2:%6.0f  epsp2:%3.2f epss2:%3.2f' %
            (cmdl.vp[1],cmdl.vs[1],cmdl.epsilonp[1],cmdl.epsilons[1]),
            xy=(anglearray[5],vs1ang[20]),xytext=(0.3,0.27),textcoords='figure fraction')
        plt.tight_layout()
        plt.savefig(pdfvve)
        if not cmdl.hideplot:
            plt.show()


    if cmdl.model != 'iso':
        ripp,rapp = Rapp(cmdl.vp[0],cmdl.vp[1],cmdl.vs[0],cmdl.vs[1],cmdl.rb[0],cmdl.rb[1],
            cmdl.delta[0],cmdl.delta[1],anglearray)

        list_ra(anglearray,rapp,ripp)
        plt.plot(anglearray,rapp,'b',label='VTI',lw=3)
        plt.plot(anglearray,ripp,'g',label= 'ISO',lw=3)
        plt.legend(loc='lower left')
        # plt.legend(loc='lower center')
        plt.title('VTI AVO')
        plt.annotate('Vp1:%6.0f Vs1:%6.0f Rb1:%3.2f d1: %3.2f' %
            (cmdl.vp[0],cmdl.vs[0],cmdl.rb[0],cmdl.delta[0]),
            xy=(anglearray[5],ripp[20]),xytext=(0.3,0.70),textcoords='figure fraction')
        plt.annotate('Vp2:%6.0f Vs2:%6.0f Rb2:%3.2f d2:%3.2f' %
            (cmdl.vp[1],cmdl.vs[1],cmdl.rb[1],cmdl.delta[1]),
            xy=(anglearray[5],ripp[20]),xytext=(0.3,0.67),textcoords='figure fraction')
        plt.ylabel('Reflection Amplitude')
        plt.xlabel('Incidence Angle in degrees')
        plt.savefig(pdfaniso0)
        if not cmdl.hideplot:
            plt.show()
        plt.plot(anglearraysqr,rapp,'b',label='VTI',lw=3)
        plt.plot(anglearraysqr,ripp,'g',label= 'ISO',lw=3)
        # plt.legend(loc='lower center')
        plt.legend(loc='lower left')
        plt.title('VTI AVO')
        plt.annotate('Vp1:%6.0f Vs1:%6.0f Rb1:%3.2f d1: %3.2f' %
            (cmdl.vp[0],cmdl.vs[0],cmdl.rb[0],cmdl.delta[0]),
            xy=(anglearraysqr[5],rapp[20]),xytext=(0.3,0.70),textcoords='figure fraction')
        plt.annotate('Vp2:%6.0f Vs2:%6.0f Rb2:%3.2f d2:%3.2f' %
            (cmdl.vp[1],cmdl.vs[1],cmdl.rb[1],cmdl.delta[1]),
            xy=(anglearraysqr[5],rapp[20]),xytext=(0.3,0.67),textcoords='figure fraction')
        plt.ylabel('Reflection Amplitude')
        plt.xlabel('Incidence Angle in radians squared')

        xi = np.linspace(anglearraysqr.min(),anglearraysqr.max())
        acfi = np.polyfit(anglearraysqr,ripp,1)
        ampii = np.polyval(acfi,xi)
        acfa = np.polyfit(anglearraysqr,rapp,1)
        ampia = np.polyval(acfa,xi)
        plt.plot(xi,ampii,c='r',lw=2,label='I_G ISO')
        plt.plot(xi,ampia,c='m',lw=2,label='I_G VTI')
        plt.annotate('ISO Intercept:%6.3f Gradient:%6.3f ' % (acfi[1],acfi[0]),
            xy=(anglearraysqr[5],rapp[20]),xytext=(0.4,0.5),textcoords='figure fraction')
        plt.annotate('VTI Intercept:%6.3f Gradient:%6.3f ' % (acfa[1],acfa[0]),
            xy=(anglearraysqr[5],rapp[20]),xytext=(0.4,0.47),textcoords='figure fraction')
        plt.legend(loc='lower left')
        # plt.legend(loc='lower center')
        plt.savefig(pdfaniso1)
        if not cmdl.hideplot:
            plt.show()

if __name__ == '__main__':
    main()
