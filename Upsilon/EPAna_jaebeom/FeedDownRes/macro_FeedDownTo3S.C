void macro_FeedDownTo3S()
{
//=========Macro generated from canvas: c3S/
//=========  (Thu Jul  7 09:30:16 2022) by ROOT version 6.24/06
   TCanvas *c3S = new TCanvas("c3S", "",0,0,700,700);
   c3S->Range(-9.25926,-12.03704,52.46914,68.20988);
   c3S->SetFillColor(0);
   c3S->SetBorderMode(0);
   c3S->SetBorderSize(2);
   c3S->SetTickx(1);
   c3S->SetTicky(1);
   c3S->SetLeftMargin(0.15);
   c3S->SetRightMargin(0.04);
   c3S->SetTopMargin(0.04);
   c3S->SetBottomMargin(0.15);
   c3S->SetFrameBorderMode(0);
   c3S->SetFrameBorderMode(0);
   
   Double_t Graph1D_y3_fx3007[2] = {
   26.5,
   34.5};
   Double_t Graph1D_y3_fy3007[2] = {
   34,
   40};
   Double_t Graph1D_y3_felx3007[2] = {
   2.5,
   5.5};
   Double_t Graph1D_y3_fely3007[2] = {
   10.63015,
   16.64332};
   Double_t Graph1D_y3_fehx3007[2] = {
   2.5,
   5.5};
   Double_t Graph1D_y3_fehy3007[2] = {
   10.63015,
   10.29563};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(2,Graph1D_y3_fx3007,Graph1D_y3_fy3007,Graph1D_y3_felx3007,Graph1D_y3_fehx3007,Graph1D_y3_fely3007,Graph1D_y3_fehy3007);
   grae->SetName("Graph1D_y3");
   grae->SetTitle(" ");
   grae->SetFillStyle(1000);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#9b43ed");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#9b43ed");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.2);
   
   TH1F *Graph_Graph1D_y33007 = new TH1F("Graph_Graph1D_y33007"," ",100,0,50);
   Graph_Graph1D_y33007->SetMinimum(0);
   Graph_Graph1D_y33007->SetMaximum(65);
   Graph_Graph1D_y33007->SetDirectory(0);
   Graph_Graph1D_y33007->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1D_y33007->SetLineColor(ci);
   Graph_Graph1D_y33007->GetXaxis()->SetTitle("p_{T}^{#varUpsilon(1S)} (GeV/c)");
   Graph_Graph1D_y33007->GetXaxis()->SetRange(1,100);
   Graph_Graph1D_y33007->GetXaxis()->CenterTitle(true);
   Graph_Graph1D_y33007->GetXaxis()->SetLabelFont(42);
   Graph_Graph1D_y33007->GetXaxis()->SetTitleSize(0.047);
   Graph_Graph1D_y33007->GetXaxis()->SetTitleOffset(1.2);
   Graph_Graph1D_y33007->GetXaxis()->SetTitleFont(42);
   Graph_Graph1D_y33007->GetYaxis()->SetTitle("Feed-down fraction (%)");
   Graph_Graph1D_y33007->GetYaxis()->CenterTitle(true);
   Graph_Graph1D_y33007->GetYaxis()->SetLabelFont(42);
   Graph_Graph1D_y33007->GetYaxis()->SetTitleSize(0.047);
   Graph_Graph1D_y33007->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph1D_y33007->GetYaxis()->SetTitleFont(42);
   Graph_Graph1D_y33007->GetZaxis()->SetLabelFont(42);
   Graph_Graph1D_y33007->GetZaxis()->SetTitleOffset(1);
   Graph_Graph1D_y33007->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph1D_y33007);
   
   grae->Draw("ap");
   
   TLegend *leg = new TLegend(0.58,0.77,0.78,0.92,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.025);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(4000);
   TLegendEntry *entry=leg->AddEntry("Graph1D_y3","LHCb 8 TeV #chi_{b}(3P) #rightarrow #varUpsilon(3S)","pe");

   ci = TColor::GetColor("#9b43ed");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#9b43ed");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   leg->Draw();
   
   TF1 *frac3pTo3s13 = new TF1("frac3pTo3s","[0]*(erf((x+[1])/[2])+1)",0,50, TF1::EAddToList::kDefault);
   frac3pTo3s13->SetFillColor(19);
   frac3pTo3s13->SetFillStyle(0);

   ci = TColor::GetColor("#990099");
   frac3pTo3s13->SetLineColor(ci);
   frac3pTo3s13->SetLineWidth(2);
   frac3pTo3s13->SetChisquare(0.109929);
   frac3pTo3s13->SetNDF(1);
   frac3pTo3s13->GetXaxis()->SetLabelFont(42);
   frac3pTo3s13->GetXaxis()->SetTitleOffset(1);
   frac3pTo3s13->GetXaxis()->SetTitleFont(42);
   frac3pTo3s13->GetYaxis()->SetLabelFont(42);
   frac3pTo3s13->GetYaxis()->SetTitleFont(42);
   frac3pTo3s13->SetParameter(0,18.18984);
   frac3pTo3s13->SetParError(0,4.23235);
   frac3pTo3s13->SetParLimits(0,0,0);
   frac3pTo3s13->SetParameter(1,-5.462886);
   frac3pTo3s13->SetParError(1,0);
   frac3pTo3s13->SetParLimits(1,-5.462886,-5.462886);
   frac3pTo3s13->SetParameter(2,11.88677);
   frac3pTo3s13->SetParError(2,0);
   frac3pTo3s13->SetParLimits(2,11.88677,11.88677);
   frac3pTo3s13->Draw("same");
   
   TPaveText *pt = new TPaveText(0.4785057,0.94,0.5214943,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText(" ");
   pt->Draw();
   c3S->Modified();
   c3S->cd();
   c3S->SetSelected(c3S);
}
