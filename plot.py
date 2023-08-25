#following command used to run plot.py where muon can be changed for bq, wboson, top_quark
# and pt can be changed to other kinematics like eta, phi maass etc.
#python plot.py ../../lhe_files/unweighted_events.root output.pdf muon pt
# how to run with command
# python plot.py ../../lhe_files/unweighted_events.root output.pdf wboson eta
# python3.11, .root is input and output.pdf will be output, wboson can be replaced with muon, top etc according 
# to the input given. eta is the variable you are plotting and can be replaced with mass, pt, phi etc.

import ROOT
import sys
if len ( sys.argv ) != 5:
    print (" USAGE : % s < input file > < output file > ") %( sys.argv [0])
    sys.exit (1)
histFileName = sys.argv [1]
plotFileName = sys.argv [2]
histo_name   = sys.argv [3]
kinematics   = sys.argv [4]
#if sys.argv[4] == "eta" :
#    kinematics == ROOT.TLatex.DrawLatex("#eta")
#elif sys.argv[4] == phi:
#    kinematics == "#phi"
print (" Reading from " , histFileName , " and writing to " , plotFileName)
histFile = ROOT.TFile.Open ( histFileName , " READ " )

# Number of histograms
histogram_names = [histo_name+f"{i}_"+kinematics for i in range(0, 4)] 

canvas = ROOT.TCanvas ( "test", "Histograms", 800, 1000 )
for idx, histogram_name in enumerate(histogram_names, start=1):
    label = histo_name+f"{idx} "+kinematics  # Construct the label
    histogram = histFile.Get(histogram_name)
    if not histogram:
        print(f"Warning: Histogram '{histogram_name}' not found in the ROOT file.")
        continue
    histogram.SetLineColor(idx)
    histogram.SetLineStyle(1)
    histogram.SetLineWidth(3)



 #   histogram = histFile.Get ( histogram_name )

#    if not histogram :
#        print (" Failed to get data histogram ")
#        sys.exit (1)

 #   histogram.SetDirectory (0)

    histogram.SetStats (0)
    histogram.SetTitle ( " " )
    histogram.SetLineWidth (2)
    histogram.GetXaxis ().SetTitle ( histo_name+f"{idx} "+kinematics+" [ GeV ] " )
    histogram.GetXaxis ().SetTitleSize (0.05)
#    histogram.GetXaxis ().SetTitleStyle (0)
    histogram.GetXaxis ().SetLabelOffset (0.005)
    histogram.GetXaxis ().SetTitleOffset (0.9)
#histogram.GetXaxis ().SetLabelSize (0.05)
    histogram.GetYaxis ().SetTitle ( " Number of events " )
    histogram.GetYaxis ().SetLabelOffset (0.005)
    histogram.GetYaxis ().SetTitleSize (0.05)
    histogram.GetYaxis ().SetTitleOffset (0.9)
    histogram.GetYaxis ().SetMaxDigits(2)

#    canvas = ROOT.TCanvas ( label, "Histograms", 1300, 1300 )
#    canvas.Clear()
    canvas.cd ()

    if idx == 1:
        histogram.Draw("h")
#    if idx == 2:
#        histogram.Draw("h")
#    if idx == 3:
#        histogram.Draw("h")
    else:
        histogram.Draw("h SAME")  # Draw on the same canvas without clearing


#    canvas.Print ( label+".png" + "[" )
#    histogram.Draw ( " h " )
#    mcHisto.Draw("same")

    legend = ROOT.TLegend (0.6 ,0.3 ,0.88 ,0.87)
    legend.AddEntry ( histogram , label, "l" )
    legend.SetLineWidth (0)
    legend.SetTextFont (1)
    legend.SetTextSize (0.05)
    legend.Draw ()

    latex = ROOT.TLatex ()
    latex.SetNDC ()
    latex.SetTextSize (0.04)
    latex.DrawLatex (0.702 ,0.905 , "#it{#sqrt{s}  = 13 TeV}" )
    latex.SetTextSize (0.04)
    latex.DrawLatex (0.19 ,0.905 , "#it{SM Four tops}" )

#    canvas.Update ()
    canvas.Print ( histo_name+f"{idx}_"+kinematics+"_"+plotFileName + "]" )

histFile.Close ()

