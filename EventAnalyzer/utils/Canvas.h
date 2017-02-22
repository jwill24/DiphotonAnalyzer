#ifndef Canvas_h
#define Canvas_h

#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <string.h>

class Canvas : public TCanvas
{
 public:
  inline Canvas(const char* name, const char* title="", bool ratio=false) :
    //TCanvas(name, "", 450, 450),
    TCanvas(name, "", 600, 600),
    fTitle(title), fTopLabel(0),
    fLeg(0), fLegX1(0.5), fLegY1(0.75),
    fRatio(ratio)
  {
    Build();
  }
  inline ~Canvas() {
    if (fLeg) delete fLeg;
    if (fTopLabel) delete fTopLabel;
  }

  inline void Prettify(TH1* obj) {
    TAxis* x = (TAxis*)obj->GetXaxis(), *y = (TAxis*)obj->GetYaxis();
    x->SetLabelFont(43); x->SetLabelSize(20);
    x->SetTitleFont(43); x->SetTitleSize(26);
    y->SetLabelFont(43); y->SetLabelSize(20);
    y->SetTitleFont(43); y->SetTitleSize(26);
    x->SetTitleColor(kBlack);
    if (fRatio) {
      x->SetTitleOffset(3.2);
      x->SetLabelOffset(0.02);
    }
    y->SetTitleOffset(1.4);

    // axis titles
    TString ttle = obj->GetTitle();
    if (ttle.Contains("\\")) {
      TObjArray* tok = ttle.Tokenize("\\");
      TString x_title = "", y_title = "", unit = "", form_spec = "";
      if (tok->GetEntries()>0) x_title = ((TObjString*)tok->At(0))->String();
      if (tok->GetEntries()>1) y_title = ((TObjString*)tok->At(1))->String();
      if (tok->GetEntries()>2) {
        unit = ((TObjString*)tok->At(2))->String();
        if (unit.Contains("?")) { // extract format specifier
          TObjArray* tok2 = unit.Tokenize("?");
          if (tok2->GetEntries()>1) {
            unit = ((TObjString*)tok2->At(0))->String();
            form_spec = ((TObjString*)tok2->At(1))->String();
          }
          else {
            unit = "";
            form_spec = ((TObjString*)tok2->At(0))->String();
          }
        }
      }
      if (!unit.IsNull() or !form_spec.IsNull()) {
        if (!unit.IsNull()) x_title = Form("%s (%s)", x_title.Data(), unit.Data());
        if (!form_spec.IsNull()) {
          TString format = Form("%%s / %%%s %%s", form_spec.Data());
          y_title = Form(format.Data(), y_title.Data(), GetBinning(obj), unit.Data());
        }
        else y_title = Form("%s / %d %s", y_title.Data(), (unsigned int)GetBinning(obj), unit.Data());
      }
      obj->GetXaxis()->SetTitle(x_title);
      obj->GetYaxis()->SetTitle(y_title);
      obj->SetTitle("");
    }
  }

  inline void RatioPlot(TH1* obj1, const TH1* obj2, const TH1* obj3, float ymin=-999., float ymax=-999.) {
    if (!fRatio) return;
    TH1* ratio1 = (TH1*)obj2->Clone(), *ratio2 = (TH1*)obj3->Clone();
    //ratio1->Sumw2(); ratio2->Sumw2();
    ratio1->Divide(obj1);
    ratio2->Divide(obj1);
    TCanvas::cd(2);
    ratio1->Draw("p");
    ratio2->Draw("p same");
    obj1->GetXaxis()->SetTitle("");
    if (ymin!=ymax) {
      ratio1->GetYaxis()->SetRangeUser(ymin, ymax);
    }
    Prettify(ratio1);
    ratio1->GetYaxis()->SetTitle("Ratios");
    TCanvas::cd();
  }

  inline void RatioPlot(TH1* obj1, const TH1* obj2, float ymin=-999., float ymax=-999.) {
    if (!fRatio) return;
    TH1* ratio = (TH1*)obj2->Clone();
    ratio->Divide(obj1);
    TCanvas::cd(2);
    ratio->Draw("p");
    obj1->GetXaxis()->SetTitle("");
    if (ymin!=ymax) {
      ratio->GetYaxis()->SetRangeUser(ymin, ymax);
    }
    Prettify(ratio);
    ratio->GetYaxis()->SetTitle("Ratio");
    TCanvas::cd();
  }

  inline void SetTopLabel(const char* lab="") {
    TCanvas::cd();
    if (strcmp(lab, "")!=0) fTitle = lab;
    if (!fTopLabel) BuildTopLabel();
    else fTopLabel->Clear();
    fTopLabel->AddText(fTitle);
    //fTopLabel->Draw();
  }

  inline TLegend* GetLegend() { return fLeg; }
  inline void SetLegendX1(double x) { fLegX1 = x; }
  inline void SetLegendY1(double y) { fLegY1 = y; }
  inline void AddLegendEntry(const TObject* obj, const char* title, Option_t* option="lpf") {
    if (!fLeg) BuildLegend();
    fLeg->AddEntry(obj, title, option);
  }

  inline void Save(const char* ext, const char* out_dir=".") {
    if (strstr(ext, "pdf")==NULL) {
      if (strstr(ext, "png")==NULL) {
        return;
      }
    }
    TCanvas::cd();
    if (fLeg and TCanvas::FindObject(fLeg)==0) fLeg->Draw();
    if (fTopLabel and TCanvas::FindObject(fTopLabel)==0) fTopLabel->Draw();
    TCanvas::SaveAs(Form("%s/%s.%s", out_dir, TCanvas::GetName(), ext));
  }

 private:
  inline void Build() {
    TCanvas::SetLeftMargin(0.14);
    TCanvas::SetTopMargin(0.06);
    TCanvas::SetRightMargin(0.1);
    TCanvas::SetBottomMargin(0.12);
    TCanvas::SetTicks(1,1);

    SetTopLabel();
    if (fRatio) DivideCanvas();
  }

  inline void DivideCanvas() {
    TCanvas::Divide(1,2);
    TPad* p1 = (TPad*)TCanvas::GetPad(1), *p2 = (TPad*)TCanvas::GetPad(2);
    p1->SetPad(0., 0.3, 1., 1.);
    p2->SetPad(0., 0.0, 1., 0.3);
    p1->SetLeftMargin(TCanvas::GetLeftMargin()); p1->SetRightMargin(TCanvas::GetRightMargin());
    p2->SetLeftMargin(TCanvas::GetLeftMargin()); p2->SetRightMargin(TCanvas::GetRightMargin());
    p1->SetTopMargin(TCanvas::GetTopMargin()+0.025); p1->SetBottomMargin(0.02);
    p2->SetTopMargin(0.02); p2->SetBottomMargin(TCanvas::GetBottomMargin()+0.225);
    p1->SetTicks(1,1); p2->SetTicks(1,1);
    p2->SetGrid(0,1);
    TCanvas::cd(1);
  }

  inline void BuildTopLabel() {
    TCanvas::cd();
    fTopLabel = new TPaveText(0.5, 0.95, 0.925, 0.96, "NB NDC");
    fTopLabel->SetFillStyle(0);
    fTopLabel->SetFillColor(0);
    fTopLabel->SetLineColor(0);
    fTopLabel->SetLineStyle(0);
    fTopLabel->SetTextFont(42);
    fTopLabel->SetTextSize(0.033);
    fTopLabel->SetTextAlign(kHAlignRight+kVAlignBottom);
  }

  inline void BuildLegend() {
    if (fLeg) return;
    if (fRatio) TCanvas::cd(1);
    fLeg = new TLegend(fLegX1, fLegY1, fLegX1+0.3, fLegY1+0.15);
    fLeg->SetFillStyle(0);
    fLeg->SetLineColor(kWhite);
    fLeg->SetLineWidth(0);
    fLeg->SetTextSize(0.04);
  }
  inline float GetBinning(const TH1* h) {
    return (h->GetXaxis()->GetXmax()-h->GetXaxis()->GetXmin())/h->GetXaxis()->GetNbins();
  }

  TString fTitle;
  TPaveText* fTopLabel;
  TLegend* fLeg;
  double fLegX1, fLegY1;
  bool fRatio;
};

#endif
