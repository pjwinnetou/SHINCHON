#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

//
// Global variables
//

TString shinchonText     = "SHINCHON";
float shinchonTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = false;
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6*2.5; 
float lumiTextOffset   = 0.01;
//float shinchonSize      = 0.75;
float shinchonTextSize      = 0.75*1.3; 
float shinchonTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

float extraOverShinchonTextSize  = 0.76;

TString lumi_pp502TeV  = "276 pb^{-1}";
TString lumi_pPb502TeV  = "34.6 nb^{-1}";
TString lumi_PbPb502TeV  = "368 #mub^{-1}";
TString lumi_PbPb502TeV_1  = "464 #mub^{-1}";
TString lumi_PbPb502TeVCent  = "368/464 #mub^{-1}";
TString lumi_sqrtS = "";

bool drawLogo      = false;

void SHINCHONLegend( TPad* pad, int iPeriod=3, int iPosX=10 ,bool isInner = true);

