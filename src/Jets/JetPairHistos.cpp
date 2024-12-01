/* **********************************************************************
 * Copyright (C) 2019-2024, Claude Pruneau, Akash Raj
 * All rights reserved.
 *
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *
 * Author: Claude Pruneau, Akash Raj,  Nov 2024
 *
 * *********************************************************************/
#include "JetPairHistos.hpp"
#include "JetHelpers.hpp"
#include "PrintHelpers.hpp"
#include "RootHelpers.hpp"
#include <math.h>
#include <cmath>

ClassImp(CAP::JetPairHistos);

namespace CAP
{

JetPairHistos::JetPairHistos()
:
HistogramGroup(),
particleDb(nullptr),
h_jet_n2_ptpt(nullptr),
h_jet_n2_phiphi(nullptr),
h_jet_n2_etaeta(nullptr),
h_jet_n2_thth(nullptr),
h_jet_n2_jtjt(nullptr),
h_jet_n2_zz(nullptr),
h_jet_n2pm_ptpt(nullptr),
h_jet_n2pm_phiphi(nullptr),
h_jet_n2pm_etaeta(nullptr),
h_jet_n2pm_thth(nullptr),
h_jet_n2pm_jtjt(nullptr),
h_jet_n2pm_zz(nullptr),
h_jet_n2mm_ptpt(nullptr),
h_jet_n2mm_phiphi(nullptr),
h_jet_n2mm_etaeta(nullptr),
h_jet_n2mm_thth(nullptr),
h_jet_n2mm_jtjt(nullptr),
h_jet_n2mm_zz(nullptr)
{
  appendClassName("JetPairHistos");
  setInstanceName("JetPairHistos");
}

JetPairHistos::JetPairHistos(const JetPairHistos & group)
:
HistogramGroup(group),
particleDb(nullptr),
h_jet_n2_ptpt(nullptr),
h_jet_n2_phiphi(nullptr),
h_jet_n2_etaeta(nullptr),
h_jet_n2_thth(nullptr),
h_jet_n2_jtjt(nullptr),
h_jet_n2_zz(nullptr),
h_jet_n2pm_ptpt(nullptr),
h_jet_n2pm_phiphi(nullptr),
h_jet_n2pm_etaeta(nullptr),
h_jet_n2pm_thth(nullptr),
h_jet_n2pm_jtjt(nullptr),
h_jet_n2pm_zz(nullptr),
h_jet_n2mm_ptpt(nullptr),
h_jet_n2mm_phiphi(nullptr),
h_jet_n2mm_etaeta(nullptr),
h_jet_n2mm_thth(nullptr),
h_jet_n2mm_jtjt(nullptr)
{
  appendClassName("JetPairHistos");
  setInstanceName("JetPairHistos");
}

JetPairHistos & JetPairHistos::operator=(const JetPairHistos & group)
{
  if (this!=&group)
  {
    HistogramGroup::operator=(group);
  }
  return *this;
}

void JetPairHistos::createHistograms()
{
  if (reportStart(__FUNCTION__)) { /* noops*/ };
  const String & bn  = getName();
  const String & ptn = getParentName();
  int nBins_pt   = configuration.getValueInt(ptn,    "nBins_pt");
  double min_pt  = configuration.getValueDouble(ptn, "min_pt");
  double max_pt  = configuration.getValueDouble(ptn, "max_pt");
  int nBins_phi  = configuration.getValueInt(ptn,    "nBins_phi");
  double min_phi = configuration.getValueDouble(ptn, "min_phi");
  double max_phi = configuration.getValueDouble(ptn, "max_phi");
  int nBins_eta  = configuration.getValueInt(ptn,    "nBins_eta");
  double min_eta = configuration.getValueDouble(ptn, "min_eta");
  double max_eta = configuration.getValueDouble(ptn, "max_eta");
  int nBins_jt  = configuration.getValueInt(ptn,    "nBins_jt");
  double min_jt = configuration.getValueDouble(ptn, "min_jt");
  double max_jt = configuration.getValueDouble(ptn, "max_jt");
  int nBins_th  = configuration.getValueInt(ptn,    "nBins_th");
  double min_th = configuration.getValueDouble(ptn, "min_th");
  double max_th = configuration.getValueDouble(ptn, "max_th");
  int nBins_z   = configuration.getValueInt(ptn,    "nBins_z");
  double min_z  = configuration.getValueDouble(ptn, "min_z");
  double max_z  = configuration.getValueDouble(ptn, "max_z");
  
  if (reportDebug(__FUNCTION__))
  {
    printCR();
    printLine();
    printValue("Parent Task Name",   ptn);
    printValue("Histo Base Name",    bn);
    printValue("nBins_pt",    nBins_pt);
    printValue("min_pt",      min_pt);
    printValue("max_pt",      max_pt);
    printValue("nBins_phi",   nBins_phi);
    printValue("min_phi",     min_phi);
    printValue("max_phi",     max_phi);
    printValue("nBins_eta",   nBins_eta);
    printValue("min_eta",     min_eta);
    printValue("max_eta",     max_eta);
    printValue("nBins_jt",    nBins_jt);
    printValue("min_jt",      min_jt);
    printValue("max_jt",      max_jt);
    printValue("nBins_th",    nBins_th);
    printValue("min_th",      min_th);
    printValue("max_th",      max_th);
    printLine();
    printCR();
  }
  h_jet_n2_ptpt    = createHistogram(createName(bn,"jet_n2_ptpt"),nBins_pt,min_pt,max_pt, nBins_pt,min_pt,max_pt,"p_{T,1}","p_{T,2}","N_{2}");
  h_jet_n2_phiphi  = createHistogram(createName(bn,"jet_n2_phiphi"),nBins_phi,min_phi,max_phi, nBins_phi,min_phi,max_phi,"#varphi_{1}","#varphi_{2}","N_{2}");
  h_jet_n2_etaeta  = createHistogram(createName(bn,"jet_n2_etaeta"),nBins_eta,min_eta,max_eta, nBins_eta,min_eta,max_eta,"#eta_{1}","#eta_{2}","N_{2}");
  h_jet_n2_thth    = createHistogram(createName(bn,"jet_n2_thth"),nBins_th,min_th,max_th, nBins_th,min_th,max_th,"#theta_{1}","#theta_{2}","N_{2}");
  h_jet_n2_jtjt    = createHistogram(createName(bn,"jet_n2_jtjt"),nBins_jt,min_jt,max_jt, nBins_jt,min_jt,max_jt,"j_{T,1}","j_{T,2}","N_{2}");
  h_jet_n2_zz      = createHistogram(createName(bn,"jet_n2_zz"),nBins_z,min_z,max_z, nBins_z,min_z,max_z,"z_{1}","z_{2}","N_{2}");
  
  h_jet_n2pp_ptpt    = createHistogram(createName(bn,"jet_n2pp_ptpt"),nBins_pt,min_pt,max_pt, nBins_pt,min_pt,max_pt,"p_{T,1}","p_{T,2}","N_{2}^{++}");
  h_jet_n2pp_phiphi  = createHistogram(createName(bn,"jet_n2pp_phiphi"),nBins_phi,min_phi,max_phi, nBins_phi,min_phi,max_phi,"#varphi_{1}","#varphi_{2}","N_{2}^{++}");
  h_jet_n2pp_etaeta  = createHistogram(createName(bn,"jet_n2pp_etaeta"),nBins_eta,min_eta,max_eta, nBins_eta,min_eta,max_eta,"#eta_{1}","#eta_{2}","N_{2}^{++}");
  h_jet_n2pp_thth    = createHistogram(createName(bn,"jet_n2pp_thth"),nBins_th,min_th,max_th, nBins_th,min_th,max_th,"#theta_{1}","#theta_{2}","N_{2}^{++}");
  h_jet_n2pp_jtjt    = createHistogram(createName(bn,"jet_n2pp_jtjt"),nBins_jt,min_jt,max_jt, nBins_jt,min_jt,max_jt,"j_{T,1}","j_{T,2}","N_{2}^{++}");
  h_jet_n2pp_zz      = createHistogram(createName(bn,"jet_n2pp_zz"),nBins_z,min_z,max_z, nBins_z,min_z,max_z,"z_{1}","z_{2}","N_{2}^{++}");

  h_jet_n2pm_ptpt    = createHistogram(createName(bn,"jet_n2pm_ptpt"),nBins_pt,min_pt,max_pt, nBins_pt,min_pt,max_pt,"p_{T,1}","p_{T,2}","N_{2}^{+-}");
  h_jet_n2pm_phiphi  = createHistogram(createName(bn,"jet_n2pm_phiphi"),nBins_phi,min_phi,max_phi, nBins_phi,min_phi,max_phi,"#varphi_{1}","#varphi_{2}","N_{2}^{+-}");
  h_jet_n2pm_etaeta  = createHistogram(createName(bn,"jet_n2pm_etaeta"),nBins_eta,min_eta,max_eta, nBins_eta,min_eta,max_eta,"#eta_{1}","#eta_{2}","N_{2}^{+-}");
  h_jet_n2pm_thth    = createHistogram(createName(bn,"jet_n2pm_thth"),nBins_th,min_th,max_th, nBins_th,min_th,max_th,"#theta_{1}","#theta_{2}","N_{2}^{+-}");
  h_jet_n2pm_jtjt    = createHistogram(createName(bn,"jet_n2pm_jtjt"),nBins_jt,min_jt,max_jt, nBins_jt,min_jt,max_jt,"j_{T,1}","j_{T,2}","N_{2}^{+-}");
  h_jet_n2pm_zz      = createHistogram(createName(bn,"jet_n2pm_zz"),nBins_z,min_z,max_z, nBins_z,min_z,max_z,"z_{1}","z_{2}","N_{2}^{+-}");

  h_jet_n2mm_ptpt    = createHistogram(createName(bn,"jet_n2mm_ptpt"),nBins_pt,min_pt,max_pt, nBins_pt,min_pt,max_pt,"p_{T,1}","p_{T,2}","N_{2}^{--}");
  h_jet_n2mm_phiphi  = createHistogram(createName(bn,"jet_n2mm_phiphi"),nBins_phi,min_phi,max_phi, nBins_phi,min_phi,max_phi,"#varphi_{1}","#varphi_{2}","N_{2}^{--}");
  h_jet_n2mm_etaeta  = createHistogram(createName(bn,"jet_n2mm_etaeta"),nBins_eta,min_eta,max_eta, nBins_eta,min_eta,max_eta,"#eta_{1}","#eta_{2}","N_{2}^{--}");
  h_jet_n2mm_thth    = createHistogram(createName(bn,"jet_n2mm_thth"),nBins_th,min_th,max_th, nBins_th,min_th,max_th,"#theta_{1}","#theta_{2}","N_{2}^{--}");
  h_jet_n2mm_jtjt    = createHistogram(createName(bn,"jet_n2mm_jtjt"),nBins_jt,min_jt,max_jt, nBins_jt,min_jt,max_jt,"j_{T,1}","j_{T,2}","N_{2}^{--}");
  h_jet_n2mm_zz      = createHistogram(createName(bn,"jet_n2mm_zz"),nBins_z,min_z,max_z, nBins_z,min_z,max_z,"z_{1}","z_{2}","N_{2}^{--}");
}

void JetPairHistos::importHistograms(TFile & inputFile)
{
  if (reportStart(__FUNCTION__)) { /* noops*/ };
  const String & bn  = getName();
  const String & ptn = getParentName();
  h_jet_n2_ptpt      = loadH2(inputFile,createName(bn,"jet_n2_ptpt"));
  h_jet_n2_phiphi    = loadH2(inputFile,createName(bn,"jet_n2_phiphi"));
  h_jet_n2_etaeta    = loadH2(inputFile,createName(bn,"jet_n2_etaeta"));
  h_jet_n2_thth      = loadH2(inputFile,createName(bn,"jet_n2_thth"));
  h_jet_n2_jtjt      = loadH2(inputFile,createName(bn,"jet_n2_jtjt"));
  h_jet_n2_zz        = loadH2(inputFile,createName(bn,"jet_n2_zz"));

  h_jet_n2pp_ptpt    = loadH2(inputFile,createName(bn,"jet_n2pp_ptpt"));
  h_jet_n2pp_phiphi  = loadH2(inputFile,createName(bn,"jet_n2pp_phiphi"));
  h_jet_n2pp_etaeta  = loadH2(inputFile,createName(bn,"jet_n2pp_etaeta"));
  h_jet_n2pp_thth    = loadH2(inputFile,createName(bn,"jet_n2pp_thth"));
  h_jet_n2pp_jtjt    = loadH2(inputFile,createName(bn,"jet_n2pp_jtjt"));
  h_jet_n2pp_zz      = loadH2(inputFile,createName(bn,"jet_n2pp_zz"));

  h_jet_n2pm_ptpt    = loadH2(inputFile,createName(bn,"jet_n2pm_ptpt"));
  h_jet_n2pm_phiphi  = loadH2(inputFile,createName(bn,"jet_n2pm_phiphi"));
  h_jet_n2pm_etaeta  = loadH2(inputFile,createName(bn,"jet_n2pm_etaeta"));
  h_jet_n2pm_thth    = loadH2(inputFile,createName(bn,"jet_n2pm_thth"));
  h_jet_n2pm_jtjt    = loadH2(inputFile,createName(bn,"jet_n2pm_jtjt"));
  h_jet_n2pm_zz      = loadH2(inputFile,createName(bn,"jet_n2pm_zz"));

  h_jet_n2mm_ptpt    = loadH2(inputFile,createName(bn,"jet_n2mm_ptpt"));
  h_jet_n2mm_phiphi  = loadH2(inputFile,createName(bn,"jet_n2mm_phiphi"));
  h_jet_n2mm_etaeta  = loadH2(inputFile,createName(bn,"jet_n2mm_etaeta"));
  h_jet_n2mm_thth    = loadH2(inputFile,createName(bn,"jet_n2mm_thth"));
  h_jet_n2mm_jtjt    = loadH2(inputFile,createName(bn,"jet_n2mm_jtjt"));
  h_jet_n2mm_zz      = loadH2(inputFile,createName(bn,"jet_n2mm_zz"));
}

void JetPairHistos::fill(PseudoJet&  jet)
{
  double jet_phi = jet.phi();
  double jet_pt  = jet.perp();
  double jet_eta = jet.pseudorapidity();
  double jet_p =   std::sqrt( ( /* px^2 */jet.px()*jet.px() ) + (/* py^2 */ jet.py()*jet.py()) + (/* pz^2 */ jet.pz()*jet.pz()) );
  /* Since we are intersed in magnitude of jT not the vector jT , i will avoid calculating jT components , This will be done down below in th following loop*/

  int pdgId1, pdgId2;
  double q1, q2;

  // Constituents of the passed Jet
  const std::vector<PseudoJet> & constituents = jet.constituents();
  for (const auto & part1 : constituents)
  {
    pdgId1 = part1.user_index();
    ParticleType * type1 = particleDb->findPdgCode(pdgId1);
    q1 = type1->getCharge();
    
    double part1_phi = part1.phi();
    double part1_pt  = part1.perp();
    double part1_eta = part1.pseudorapidity();
    double part1_z   = part1_pt/jet_pt;

    double part1_p = std::sqrt( (/*px^2*/ part1.px()*part1.px()) + (/*py^2*/ part1.py()*part1.py()) +(/*pz^2*/ part1.pz()*part1.pz()) );
    double part1_jt  = std::sqrt( (/* particle_p^2*/part1_p*part1_p) - (/* jet_p^2*/jet_p*jet_p) );    /* using loging jT^2 = p_particle^2 + p_Jet^2*/ 
    double part1_th  = std::acos(jet_p/part1_p); /* usnig arccos(Jet's_p/particle_p) to finf the theta between the constituent particle vector from the Jets axix*/

    //calculateJtTheta(jet_px,jet_py,jet_pz,jet_p,part1_px,part1_py,part1_pz,part1_p,part1_jt,part1_th);
    for (const auto & part2 : constituents)
    {
      if (part1==part2) continue;
      
      pdgId2 = part2.user_index();
      ParticleType * type2 = particleDb->findPdgCode(pdgId2);
      q2 = type2->getCharge();

      double part2_phi = part2.phi();
      double part2_pt  = part2.perp();
      double part2_eta = part2.pseudorapidity();
      double part2_z   = part2_pt/jet_pt;

      double part2_p = std::sqrt( (/*px^2*/ part2.px()*part2.px()) + (/*py^2*/ part2.py()*part2.py()) +(/*pz^2*/ part2.pz()*part2.pz()) );
      double part2_jt  = std::sqrt( (/* particle_p^2*/part2_p*part2_p) - (/* jet_p^2*/jet_p*jet_p) );    /* using loging jT^2 = p_particle^2 + p_Jet^2*/ 
      double part2_th  = std::acos(jet_p/part2_p); /* usnig arccos(Jet's_p/particle_p) to finf the theta between the constituent particle vector from the Jets axix*/

      //calculateJtTheta(jet_px,jet_py,jet_pz,jet_p,part1_px,part2_py,part2_pz,part2_p,part2_jt,part2_th);
      
      h_jet_n2_ptpt->Fill(part1_pt,part2_pt);
      h_jet_n2_phiphi->Fill(part1_phi,part2_phi);
      h_jet_n2_etaeta->Fill(part1_eta,part2_eta);
      h_jet_n2_thth->Fill(part1_th,part2_th);
      h_jet_n2_jtjt->Fill(part1_jt,part2_jt);
      h_jet_n2_zz->Fill(part1_z,part2_z);
      
      if (q1>0 && q2>0)
      {
        h_jet_n2pp_ptpt->Fill(part1_pt,part2_pt);
        h_jet_n2pp_phiphi->Fill(part1_phi,part2_phi);
        h_jet_n2pp_etaeta->Fill(part1_eta,part2_eta);
        h_jet_n2pp_thth->Fill(part1_th,part2_th);
        h_jet_n2pp_jtjt->Fill(part1_jt,part2_jt);
        h_jet_n2pp_zz->Fill(part1_z,part2_z);
      }
      else if (q1>0 && q2<0)
      {
        h_jet_n2pm_ptpt->Fill(part1_pt,part2_pt);
        h_jet_n2pm_phiphi->Fill(part1_phi,part2_phi);
        h_jet_n2pm_etaeta->Fill(part1_eta,part2_eta);
        h_jet_n2pm_thth->Fill(part1_th,part2_th);
        h_jet_n2pm_jtjt->Fill(part1_jt,part2_jt);
        h_jet_n2pm_zz->Fill(part1_z,part2_z);
      }
      else if (q1<0 && q2>0)
      {
        h_jet_n2pm_ptpt->Fill(part2_pt,part1_pt);
        h_jet_n2pm_phiphi->Fill(part2_phi,part1_phi);
        h_jet_n2pm_etaeta->Fill(part2_eta,part1_eta);
        h_jet_n2pm_thth->Fill(part2_th,part1_th);
        h_jet_n2pm_jtjt->Fill(part2_jt,part1_jt);
        h_jet_n2pm_zz->Fill(part2_z,part1_z);
      }
      else if (q1<0 && q2<0)
      {
        h_jet_n2mm_ptpt->Fill(part1_pt,part2_pt);
        h_jet_n2mm_phiphi->Fill(part1_phi,part2_phi);
        h_jet_n2mm_etaeta->Fill(part1_eta,part2_eta);
        h_jet_n2mm_thth->Fill(part1_th,part2_th);
        h_jet_n2mm_jtjt->Fill(part1_jt,part2_jt);
        h_jet_n2mm_zz->Fill(part1_z,part2_z);
      }
    }
  }
}

void JetPairHistos::scaleHistograms(double scale)
{
  h_jet_n2_ptpt->Scale(scale);
  h_jet_n2_phiphi->Scale(scale);
  h_jet_n2_etaeta->Scale(scale);
  h_jet_n2_thth->Scale(scale);
  h_jet_n2_jtjt->Scale(scale);
  h_jet_n2_zz->Scale(scale);
  
  h_jet_n2pp_ptpt->Scale(scale);
  h_jet_n2pp_phiphi->Scale(scale);
  h_jet_n2pp_etaeta->Scale(scale);
  h_jet_n2pp_thth->Scale(scale);
  h_jet_n2pp_jtjt->Scale(scale);
  h_jet_n2pp_zz->Scale(scale);
  
  h_jet_n2pm_ptpt->Scale(scale);
  h_jet_n2pm_phiphi->Scale(scale);
  h_jet_n2pm_etaeta->Scale(scale);
  h_jet_n2pm_thth->Scale(scale);
  h_jet_n2pm_jtjt->Scale(scale);
  h_jet_n2pm_zz->Scale(scale);

  h_jet_n2mm_ptpt->Scale(scale);
  h_jet_n2mm_phiphi->Scale(scale);
  h_jet_n2mm_etaeta->Scale(scale);
  h_jet_n2mm_thth->Scale(scale);
  h_jet_n2mm_jtjt->Scale(scale);
  h_jet_n2mm_zz->Scale(scale);
  
}

} // namespace CAP
