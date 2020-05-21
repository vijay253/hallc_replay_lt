// Reads the points from a file and produces a simple graph.

int macro2()
{
  auto c = new TCanvas();c->SetGrid();

  TGraphErrors graph_expected("input.txt", "%lg %lg %lg");
                               
  graph_expected.SetTitle( "Consistency of calibration constants for PMT1;" "Run Number;" " Calibration constant");			  			 			 
  graph_expected.SetFillColor(kYellow);
  graph_expected.DrawClone("E3AL"); // E3 draws the band

  TGraphErrors graph("input.txt","%lg %lg %lg");
  graph.SetMarkerStyle(kCircle);
  graph.SetFillColor(0);
  graph.DrawClone("PESame");

  // Draw the Legend
  TLegend leg(.1,.7,.3,.9,"Lab. Lesson 2");
  leg.SetFillColor(0);
  leg.AddEntry(&graph_expected,"First Guess");
  leg.AddEntry(&graph,"Second Guess");
  leg.DrawClone("Same");
  graph.Print();
  return 0;
}
