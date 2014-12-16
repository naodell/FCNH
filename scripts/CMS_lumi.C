#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

TString cmsText     = "CMS";
float cmsTextFont   = 42;

bool writeExtraText = false;
TString extraText   = "Preliminary";
float extraTextFont = 52;

TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";

bool outOfFrame    = false;
bool drawLogo      = false;

void CMS_lumi( TPad* pad, int iPeriod=3, 
	       int iPosX = 0, 
	       float alpha=0.15, float beta=0.05, 
	       float delta=0.90 );

void 
CMS_lumi( TPad* pad, int iPeriod, 	  
	  int iPosX, 
	  float alpha, float beta, float delta )
{            
  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  float e = 0.025;

  pad->cd();

  TString lumiText;
  if( iPeriod==1 )
    {
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==2 )
    {
      lumiText += lumi_8TeV;
      lumiText += " (8 TeV)";
    }
  else if( iPeriod==3 ) 
    {
      lumiText = lumi_8TeV; 
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==4 )
    {
      lumiText += lumi_13TeV;
      lumiText += " (13 TeV)";
    }
  else if ( iPeriod==7 )
    { 
      if( outOfFrame ) lumiText += "#scale[0.85]{";
      lumiText += lumi_13TeV; 
      lumiText += " (13 TeV)";
      lumiText += " + ";
      lumiText += lumi_8TeV; 
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
      if( outOfFrame) lumiText += "}";
    }
   
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float lumiTextSize   = 0.6*t;
  float lumiTextOffset = 0.2*t;
  float cmsTextSize    = 0.99*t;
  float cmsTextOffset  = 0.1*t;
  float extraTextSize  = 1.0*lumiTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset,lumiText);

  if( outOfFrame )
    {
      latex.SetTextFont(42);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize);    
      latex.DrawLatex(l,1-t+cmsTextOffset,cmsText);
    }
  
  pad->cd();

  float posX_;
  if( iPosX==0 )
    {
      posX_ =   l + alpha*(1-l-r);
    }
  else if( iPosX==1 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX==2 )
    {
      posX_ =  1-r - alpha*(1-l-r);
    }
  float posY_ = 1-t - beta*(1-t-b);
  if( !outOfFrame )
    {
      if( drawLogo )
	{
	  posX_ =   l + 0.045*(1-l-r)*W/H;
	  posY_ = 1-t - 0.045*(1-t-b);
	  float xl_0 = posX_;
	  float yl_0 = posY_ - 0.15;
	  float xl_1 = posX_ + 0.15*H/W;
	  float yl_1 = posY_;
	  TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
	  TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
	  pad_logo->Draw();
	  pad_logo->cd();
	  CMS_logo->Draw("X");
	  pad_logo->Modified();
	  pad->cd();
	}
      else
	{
	  latex.SetTextFont( cmsTextFont );
	  latex.SetTextSize( cmsTextSize );
	  latex.SetTextAlign( 23 );
	  latex.DrawLatex( posX_, posY_, cmsText );
	  if( writeExtraText ) 
	    {
	      latex.SetTextFont( extraTextFont );
	      latex.SetTextAlign( 23 );
	      latex.SetTextSize( extraTextSize );
	      latex.DrawLatex( posX_, posY_- delta*cmsTextSize, extraText );
	    }
	}
    }
  else if( writeExtraText )
    {
      latex.SetTextFont( extraTextFont );
      latex.SetTextSize( extraTextSize );
      latex.SetTextAlign( 23 );
      latex.DrawLatex( posX_, posY_, extraText );      
    }
  return;
}
