#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"

void plot(){

// The values and the errors on the Y axis

const int n_points=10;
double x_vals[n_points]= {1,2,3,4,5,6,7,8,9,10};
double y_vals[n_points]= {6,12,14,20,22,24,35,45,44,53};
double y_errs[n_points]= {5,5,4.7,4.5,4.2,5.1,2.9,4.1,4.8,5.43};
double y_vals1[n_points]= {2,4,10,18,20,21,30,40,41,50};

// Instance of the graph

TGraphErrors graph(n_points,x_vals,y_vals,nullptr,y_errs);
graph.SetTitle("Comparison of colibration constants; Run No ; Calibration constants");

TGraphErrors graph1(n_points,x_vals,y_vals1,nullptr,y_errs);


// Make the plot estetically better

graph.SetMarkerStyle(kOpenCircle);
graph.SetMarkerColor(kBlue);
graph.SetLineColor(kBlue);

// The canvas on which we'll draw the graph

auto mycanvas =new TCanvas();

// Draw the graph !

graph.DrawClone("APE");
graph1.DrawClone("same");
 mycanvas->SaveAs("plot.pdf");

// Define a linear function

TF1 f("Linear law","[0]+x*[1]",.5,10.5);
// Let's make the function line nicer

f.SetLineColor(kRed); f.SetLineStyle(2);
// Fit it to the graph and draw it
graph.Fit(&f);
f.DrawClone("Same");
// Build and Draw a legend
TLegend leg(.1,.7,.3,.9,"Lab. Lesson 1");
leg.SetFillColor(0);
graph.SetFillColor(0);
leg.AddEntry(&graph,"Exp. Points");
leg.AddEntry(&f,"Th. Law");
leg.DrawClone("Same");
// Draw an arrow on the canvas
TArrow arrow(8,8,6.2,23,0.02,"|>");
arrow.SetLineWidth(2);
arrow.DrawClone();
// Add some text to the plot
TLatex text(8.2,7.5,"#splitline{Maximum}{Deviation}");
text.DrawClone();
mycanvas->Print("graph_with_law.pdf");
}
int main(){
  plot();}
