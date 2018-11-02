# keep ROOT TApplication from grabbing -h flag
from ROOT import PyConfig
PyConfig.IgnoreCommandLineOptions = True
from ROOT import *

from pdeTest import *
import sys,os
import csv
import argparse
import numpy




if len(sys.argv)<2:
    print "No input file/template given"
    print "Usage: pdeFitParams.py filename or template"
    sys.exit()

parser = argparse.ArgumentParser(description='pdeFitParams analyzer') 
parser.add_argument('files', nargs='*', help="Give file paths to specify files.")
parser.add_argument("-p", "--plotAll", default=None, help="show the plots",
                    action="store_true")
parser.add_argument('-c', '--calib', type=float, default=-1,
                        help="Ref. current to photons calibration")

args = parser.parse_args()
calib=args.calib

vbias=[]
refcurrent=[]
picorange=[]
npeval=[]
vex=[]
v_ov=[]
params=[]
freq=10000
ratio = 4690
resp = .427085
wlen = 640E-9
h = 6.62607E-34
c = 2.99E8
phot_e = h*c/wlen
nphot_ref=[]
nphot_test=[]
pde=[]
gain_arr=[]
noise_arr=[]
enf_arr=[]



print "Ref. current to photons calibration",calib
tgGain=TGraph()
tgGain.SetTitle("Gain vs. Vbias;Vbias [V];Gain [arbitrary units]")

tgPDE=TGraph()
tgPDE.SetTitle("PDE vs V_ex; V_ex [V]; PDE [%]")

tgVbr=TGraph()

tgVex=TGraph()

tgNoise=TGraph()

tgENF=TGraph()

for fn in args.files:
    print fn
    parsed=os.path.basename(fn).replace(".root","").split("-")
    print parsed
    dac = "{:.0f}".format((float(parsed[5].replace("dac",""))))
    vbias.append(float(parsed[6].replace("_",".").replace("V","")))
    if "pa" in parsed[7]:
        refcurrent.append(float(parsed[7].replace("_",".").replace("pa",""))/1000)
    elif:
        refcurrent.append(float(parsed[7].replace("_",".").replace("na","")))
    else:
        refcurrent.append(float(parsed[7].replace("_",".").replace("Na","")))
    chip_id = parsed[3]
    chip_ch = parsed[4]
print vbias, "V"
print refcurrent, "nA"



csvfile=open('parameters.csv','w')
writer=csv.writer(csvfile)


for i in range(len(args.files)):

    tf=TFile(args.files[i])
    hLight=tf.Get("hpulses1")
    hDark=tf.Get("hpulses0")
    hRange=tf.Get("hRange")
    picorange.append(hRange.GetBinContent(1))

    if (args.plotAll):
        hLight.Draw()
        hDark.Draw("same")

    ana=PhDAnalyzier(hLight.Clone(),hDark.Clone())
    npe=ana.CalcNpe()
        
    ana.FitPhD() # do a nice fit to the peaks

    if (args.plotAll):
        screenY=TGClient.Instance().GetDisplayHeight()
        c1=TCanvas("results","results",int(screenY*.75),int(screenY*.75))
        c1.Divide(1,2)
        c1.cd(1)
        ana.hPhD.Draw()
        c1.cd(2)
        ana.hPhD0.Draw()

    print "Mean NPE detectected",npe
    
    if (args.plotAll): raw_input("Press Enter to continue...")

    gain=ana.GetGain()*picorange[i]/100  # scale all gains to 100 mV range
    tgGain.SetPoint(i,vbias[i],gain)
        
    row=[sys.argv[1],npe,ana.GetNoise(),ana.GetGain(),ana.GetENF(),tf.Get("hRange").GetBinContent(1)]
    writer.writerow(row)

    npeval.append(npe) #Fill NPE array
    gain_arr.append(gain) #Fill Gain array

    noise_arr.append(ana.GetNoise())
    enf_arr.append(ana.GetENF())
    
    
print "gain", gain_arr   

csvfile.close()

tcan=TCanvas("fitparams","Fit params")
tgGain.Fit("pol1")
tgGain.Draw("ALP*")

func = tgGain.GetFunction("pol1")

#Determine x-intercept -> V_breakdown
p0 = func.GetParameter(0);
p1 = func.GetParameter(1);
vbr =  -round(p0/p1,2)
ylast = round(p1*(vbias[len(vbias)-1])+p0,2)

#Define points of extrapolated line
vbr_fit = [(vbr,0), vbias[len(vbias)-1],ylast]

vbr_x =[vbr, vbias[(len(vbias)-1)]]
vbr_y =[0, ylast]

print "vbry", vbr_y
print "vbrx", vbr_x

#Set points to draw extrapolated V_br line
for i in range(len(vbr_x)):
    tgVbr.SetPoint(i,vbr_x[i],vbr_y[i])
    
tgVbr.Draw("AlP")

tcan.Draw()

#Gain vs. Voltage graph settings
tgGain.GetXaxis().SetLimits(vbr*.99,vbias[len(vbias)-1])
tgGain.GetYaxis().SetRangeUser(0,1.1*max(gain_arr))
tgGain.Draw("ALP*")
tgVbr.GetXaxis().SetLimits(vbr*.99,vbias[len(vbias)-1])
tgVbr.GetYaxis().SetRangeUser(0,1.1*max(gain_arr))
tgVbr.Draw("same")
tgVbr.SetLineColor(2)
tgVbr.SetLineWidth(2)

text = TLatex()
text.SetTextAlign(13)
text.DrawLatex(vbr,-500,Form("V_br = %g V" % vbr))

text.DrawClone()
tcan.Modified()

#Determine number of photons at test/ref per each ref. curr.
for i in range(len(args.files)):
    vex.append(vbias[i] - vbr)
    nphot_ref.append((refcurrent[i]*1E-9)/(resp*freq*phot_e))
    nphot_test.append(nphot_ref[i]/ratio)

#From number of photons determine NPE at DUT
##Then plot various parameters
for i in range(len(args.files)):
    pde.append(npeval[i]/nphot_test[i])
    #set points for plotting, PDE, Gain,ENF,  and Noise v V_ex
    tgPDE.SetPoint(i,vex[i], pde[i]*100)
    tgVex.SetPoint(i,vex[i], gain_arr[i])
    tgNoise.SetPoint(i,vex[i],noise_arr[i])
    tgENF.SetPoint(i,vex[i],enf_arr[i])

tcanpde=TCanvas("pde", "PDE vs V_ex")
#tgPDE.Fit("pol3")
tgPDE.Draw("ALP*")
tcanpde.Draw()

tcanvex=TCanvas("vex", "Gain vs V_ex")
tgVex.Draw("ALP*")
tcanvex.Draw()

tgVex.GetXaxis().SetTitle("V_ex [V]")
tgVex.GetYaxis().SetTitle("Gain [A. U.]")
tgVex.SetTitle("Gain vs V_ex")

##Putting all plots onto a single split canvas


allcan = TCanvas("allcan", '{0} {1} PDE Plots ({2} dac)'.format(chip_id, chip_ch, dac),2000,1700)

#Canvas Title
allcanLabel = TPaveLabel(.1,.96,.9,.99, '{0} {1} PDE Plots ({2} dac)'.format(chip_id, chip_ch,dac))
allcanLabel.Draw()
graphPad = TPad("Graphs", "Graphs",0.06,0.05,0.95,.95)
graphPad.Draw()
graphPad.cd()
graphPad.Divide(3,2)


#plot pde vs V_ex
graphPad.cd(1)
tgPDE.Draw()

#plot gain v V_ex
graphPad.cd(2)
tgVex.Draw()

#plot Gain v V_bias
labely = -0.075*max(gain_arr)
graphPad.cd(3)
tgGain.Draw()
tgVbr.Draw("same")
text = TLatex()
text.SetTextAlign(13)
text.DrawLatex(vbr,labely,Form("V_br = %g V" % vbr))

#Plot Noise v V_ex
graphPad.cd(4)
tgNoise.Draw("ALP*")
tgNoise.GetXaxis().SetTitle("V_ex [V]")
tgNoise.GetYaxis().SetTitle("Noise [A. U.]")
tgNoise.SetTitle("Noise vs V_ex")

#Plot ENF vs V_ex
graphPad.cd(5)
tgENF.Draw("ALP*")
tgENF.GetXaxis().SetTitle("V_ex [V]")
tgENF.GetYaxis().SetTitle("ENF [A. U]")
tgENF.SetTitle("ENF vs V_ex")
#graphPad.cd(5).SetLeftMargin(.15)

for i in range(1,6):
    graphPad.cd(i).SetLeftMargin(.15)
                 
print "vbr is: ", vbr
print "dac is: ", dac




raw_input("Press Enter to continue...")
