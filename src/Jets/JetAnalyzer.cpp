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
#include "JetAnalyzer.hpp"
#include "JetHistos.hpp"
#include "JetSingleHistos.hpp"
#include "JetPairHistos.hpp"
#include "PrintHelpers.hpp"
#include "fastjet/ClusterSequence.hh"

ClassImp(CAP::JetAnalyzer);

namespace CAP
{

JetAnalyzer::JetAnalyzer()
:
EventTask()
{
  appendClassName("JetAnalyzer");
  setInstanceName("JetAnalyzer");
  setName("JetAnalyzer");
  setTitle("JetAnalyzer");
}

JetAnalyzer::JetAnalyzer(const JetAnalyzer & task)
:
EventTask(task)
{   }

JetAnalyzer & JetAnalyzer::operator=(const JetAnalyzer & task)
{
  if (this!=&task)
    {
    EventTask::operator=(task);
    }
  return *this;
}


void JetAnalyzer::setDefaultConfiguration()
{
  EventTask::setDefaultConfiguration();
  addProperty("HistogramBaseName","Jet");
  addProperty("JetRadius",      0.4);
  addProperty("JetPtMin",       10.0);
  addProperty("nBins_jet_n1",   100);
  addProperty("min_jet_n1",     0.0);
  addProperty("max_jet_n1",     100.0);
  addProperty("nBins_jet_p",    100);
  addProperty("min_jet_p",      0.0);
  addProperty("max_jet_p",      100.0);
  addProperty("nBins_jet_pt",   100);
  addProperty("min_jet_pt",     0.0);
  addProperty("max_jet_pt",     100.0);
  addProperty("nBins_jet_phi",  32);
  addProperty("min_jet_phi",    0.0);
  addProperty("max_jet_phi",    CAP::Math::twoPi());
  addProperty("nBins_jet_eta",   20);
  addProperty("min_jet_eta",     -2.0);
  addProperty("max_jet_eta",     2.0);
  addProperty("nBins_jet_netQ",  20);
  addProperty("min_jet_netQ",    -10.0);
  addProperty("max_jet_netQ",    10.0);
  addProperty("nBins_p",    100);
  addProperty("min_p",      0.0);
  addProperty("max_p",      10.0);
  addProperty("nBins_pt",   100);
  addProperty("min_pt",     0.0);
  addProperty("max_pt",     10.0);
  addProperty("nBins_phi",  36);
  addProperty("min_phi",    0.0);
  addProperty("max_phi",    CAP::Math::twoPi());
  addProperty("nBins_eta",  40);
  addProperty("min_eta",    -2.0);
  addProperty("max_eta",    2.0);
  addProperty("nBins_jt",   50);
  addProperty("min_jt",     0.0);
  addProperty("max_jt",     5.0);
  addProperty("nBins_th",   20);
  addProperty("min_th",     0.0);
  addProperty("max_th",     CAP::Math::pi()/4.0);
  addProperty("nBins_z",   50);
  addProperty("min_z",     0.0);
  addProperty("max_z",     1.0);
}

void JetAnalyzer::configure()
{
  EventTask::configure();
}

void JetAnalyzer::initialize()
{
  if (reportStart(__FUNCTION__)) {/* no ops */}
  EventTask::initialize();

  particleDb = Manager<ParticleDb>::getObjectAt(0);
 
  clearSets();
  addSet("JetHistos");
  addSet("JetSingleHistos");
  addSet("JetPairHistos");
  createHistograms();
  if (reportEnd(__FUNCTION__)) {/* no ops */}
}

void JetAnalyzer::createHistograms()
{
  if (reportStart(__FUNCTION__)) {/* no ops */};
  String bn = getValueString("HistogramBaseName");
  jetRadius = getValueDouble("JetRadius");
  jetPtMin  = getValueDouble("JetPtMin");
  JetHistos * jetHistos;
  JetSingleHistos * jetSingleHistos;
  JetPairHistos * jetPairHistos;

  for (auto & eventFilter : eventFilters)
    {
    String efn = eventFilter->getName();
    for (auto & jetFilter : jetFilters)
      {
      String jfn = jetFilter->getName();

      jetHistos = new JetHistos();
      jetHistos->setName(createName(bn,efn,jfn));
      jetHistos->setConfiguration(configuration);
      jetHistos->setParentTask(this);
      jetHistos->createHistograms();
      jetHistos->setParticleDb(getParticleDb());
      addGroupInSet(0,jetHistos);

      jetSingleHistos = new JetSingleHistos();
      jetSingleHistos->setName(createName(bn,efn,jfn));
      jetSingleHistos->setConfiguration(configuration);
      jetSingleHistos->setParentTask(this);
      jetSingleHistos->createHistograms();
      jetSingleHistos->setParticleDb(getParticleDb());
      addGroupInSet(1,jetSingleHistos);

      jetPairHistos = new JetPairHistos();
      jetPairHistos->setName(createName(bn,efn,jfn));
      jetPairHistos->setConfiguration(configuration);
      jetPairHistos->setParentTask(this);
      jetPairHistos->createHistograms();
      jetPairHistos->setParticleDb(getParticleDb());
      addGroupInSet(2,jetPairHistos);
      }
    }
  if (reportEnd(__FUNCTION__)) {/* no ops */};
}

void JetAnalyzer::execute()
{
  Event & event = *Manager<Event>::getObjectAt(0);
  std::vector<Particle*> & particles = event.getParticles();
  if (!analyzeThisEvent(event,eventFilters,eventFilterAccepted)) return;
  std::vector<PseudoJet> pseudoJetsInput;
  pseudoJetsInput.clear();

//  printCR;
//  printLine();
//  printString("New Event");
//  printValue("particles.size()",particles.size());
  //int k=0;
  for (auto & particle : particles)
    {
    int pid = particle->getType().getPdgCode();
    LorentzVector & momentum = particle->getMomentum();
//    printValue("k",k);
//    printValue("px",momentum.Px());
//    printValue("py",momentum.Py());
//    printValue("pz",momentum.Pz());
//    printValue("e",momentum.E());
//    printValue("pid",pid);
    PseudoJet pseudoJet (momentum.Px(),momentum.Py(),momentum.Pz(),momentum.E());
    pseudoJet.set_user_index(pid);
    if(pseudoJet.eta()>=-3 && pseudoJet.eta()<=3){
      pseudoJetsInput.push_back(pseudoJet);
    }else{
      continue;
    }
    //k++;
    }
//  printValue("jetRadius",jetRadius);
//  printValue("jetPtMin",jetPtMin);

  JetDefinition jetDef(antikt_algorithm,jetRadius);
  ClusterSequence * clusterSequence = new ClusterSequence(pseudoJetsInput, jetDef);
  std::vector<PseudoJet> clusteredJets = sorted_by_pt( clusterSequence->inclusive_jets(jetPtMin) );
//printJets (clusteredJets);
//  printValue("jetRadius",jetRadius);
//  printValue("clusteredJets.size()",clusteredJets.size());


  unsigned int iEventFilter = 0;
  for (auto accepted : eventFilterAccepted)
    {
    if (!accepted)
      {
      printString("EVENT NOT ACCEPTED");
      continue;
      }
    //printString("EVENT ACCEPTED");


    unsigned int  baseSingle   = iEventFilter*nJetFilters;
    for (unsigned int iJetFilter=0; iJetFilter<nJetFilters;iJetFilter++)
      {
      int index = baseSingle+iJetFilter;
      JetHistos & jetHistos = (JetHistos &) getGroupAt(0,index);
      JetSingleHistos & jetSingleHistos = (JetSingleHistos &) getGroupAt(1,index);
      JetPairHistos & jetPairHistos     = (JetPairHistos &) getGroupAt(2,index);
      for (auto jet : clusteredJets)
        {
        if (jetFilters[iJetFilter]->accept(jet))
          {
          jetHistos.fill(jet);
          jetSingleHistos.fill(jet);
          jetPairHistos.fill(jet);
          }
        }
      }
    iEventFilter++;
    }
  TaskAccountant::increment();
  delete clusterSequence;
}

//!
//!Jet histograms are scaled by the number of jets accepted - instead of the number of events accepted.
//!
void JetAnalyzer::scaleHistograms()
{
  if (reportStart(__FUNCTION__)) { /* no ops */ };

  int nEventFilters = eventFilters.size();
  int nJetFilters   = jetFilters.size();
  for (int iEventFilter=0; iEventFilter<nEventFilters; iEventFilter++)
  {
    unsigned int  baseSingle   = iEventFilter*nJetFilters;
    for (int iJetFilter=0; iJetFilter<nJetFilters; iJetFilter++)
    {
      int index = baseSingle+iJetFilter;
      JetHistos & jetHistos = (JetHistos &) getGroupAt(0,index);
      double jetAccepted = jetHistos.getAcceptedJets();
      double scale = 1.0;
      if (jetAccepted>1)
        scale = 1.0/jetAccepted;
      jetHistos.scaleHistograms(scale);
      JetSingleHistos & jetSingleHistos = (JetSingleHistos &) getGroupAt(1,index);
      jetSingleHistos.scaleHistograms(scale);
      JetPairHistos & jetPairHistos     = (JetPairHistos &) getGroupAt(2,index);
      jetPairHistos.scaleHistograms(scale);
    }
  }
  if (reportEnd(__FUNCTION__)) {/* no ops */};
}

void JetAnalyzer::printJets (const vector<fastjet::PseudoJet> & jets)
{
  // sort jets into increasing pt
  vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);

  // label the columns
  printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity",
         "phi", "pt", "n constituents");

  // print out the details for each jet
  for (unsigned int i = 0; i < sorted_jets.size(); i++)
    {
    // the following is not super efficient since it creates an
    // intermediate constituents vector
    int n_constituents = sorted_jets[i].constituents().size();
    printf("%5u %15.8f %15.8f %15.8f %8u\n",
           i, sorted_jets[i].rap(), sorted_jets[i].phi(),
           sorted_jets[i].perp(), n_constituents);
  }
}

} // namespace CAP
