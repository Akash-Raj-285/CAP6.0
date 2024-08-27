/* **********************************************************************
 * Copyright (C) 2019-2022, Claude Pruneau, Victor Gonzalez, Sumit Basu
 * All rights reserved.
 *
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *
 * Author: Claude Pruneau,   04/01/2022
 *
 * *********************************************************************/
#include "PtPtHistos.hpp"
#include "PrintHelpers.hpp"
#include "RootHelpers.hpp"

ClassImp(CAP::PtPtHistos);

namespace CAP
{

PtPtHistos::PtPtHistos()
:
HistogramGroup(),
multName("MULT"),
h_mult(),
h_f1_vsMult(),
h_f2_vsMult(),
h_f3_vsMult(),
h_f4_vsMult(),
h_Q1_vsMult(),
h_Q1Q1_vsMult(),
h_Q1Q1Q1_vsMult(),
h_Q1Q1Q1Q1_vsMult(),
h_Q2_vsMult(),
h_Q2Q1_vsMult(),
h_Q2Q1Q1_vsMult(),
h_Q2Q2_vsMult(),
h_Q3_vsMult(),
h_Q3Q1_vsMult(),
h_Q4_vsMult(),
particleFilters()
{
  appendClassName("PtPtHistos");
  setInstanceName("PtPt");
}

PtPtHistos::PtPtHistos(const PtPtHistos & group)
:
HistogramGroup(group),
multName(group.multName),
h_mult(),
h_f1_vsMult(),
h_f2_vsMult(),
h_f3_vsMult(),
h_f4_vsMult(),
h_Q1_vsMult(),
h_Q1Q1_vsMult(),
h_Q1Q1Q1_vsMult(),
h_Q1Q1Q1Q1_vsMult(),
h_Q2_vsMult(),
h_Q2Q1_vsMult(),
h_Q2Q1Q1_vsMult(),
h_Q2Q2_vsMult(),
h_Q3_vsMult(),
h_Q3Q1_vsMult(),
h_Q4_vsMult(),
particleFilters()
{
  cloneAll(group);
}

PtPtHistos & PtPtHistos::operator=(const PtPtHistos & group)
{
  if (this!=&group)
    {
    HistogramGroup::operator=(group),
    cloneAll(group);
    }
  return *this;
}

void PtPtHistos::cloneAll(const PtPtHistos & group  __attribute__((unused)))
{
  if (reportStart(__FUNCTION__))  { /* no ops */ }
  if (reportEnd(__FUNCTION__))  { /* no ops */ }
}


// for now use the same boundaries for eta and y histogram
void PtPtHistos::createHistograms()
{
  if (reportStart(__FUNCTION__))  { /* no ops */ }
  const String & bn  = getName();
  const String & ptn = getParentName();
  const String & ppn = getParentPathName();
  multName = configuration.getValueInt(ppn,"MultName");
  int nBins_mult   = configuration.getValueInt(ppn,"nBins_mult");
  double min_mult  = configuration.getValueInt(ppn,"Min_mult");
  double max_mult  = configuration.getValueInt(ppn,"Max_mult");

  if (reportInfo(__FUNCTION__))
    {
    printCR();
    printValue("Parent Task Name",ptn);
    printValue("Parent Path Name",ppn);
    printValue("Histo Base Name.",bn);
    printValue("multName",multName);
    printValue("nBins_mult",nBins_mult);
    printValue("min_mult",min_mult);
    printValue("max_mult",max_mult);
    }
  String xTitle = "M";
  h_mult = createHistogram(createName(bn,"mult"),nBins_mult,min_mult,max_mult,xTitle,"Counts");
  // use multiple filters -- but diagonnally only
  // no mixed species for now.
  for (auto & particleFilter1 : particleFilters)
    {
    String pfn = particleFilter1->getName();
    h_f1_vsMult.push_back( createProfile(createName(bn,pfn,"n1"),nBins_mult,min_mult,max_mult,xTitle,"n_{1}") );
    h_f2_vsMult.push_back( createProfile(createName(bn,pfn,"n2"),nBins_mult,min_mult,max_mult,xTitle,"n_{2}") );
    h_f3_vsMult.push_back( createProfile(createName(bn,pfn,"n3"),nBins_mult,min_mult,max_mult,xTitle,"n_{3}") );
    h_f4_vsMult.push_back( createProfile(createName(bn,pfn,"n4"),nBins_mult,min_mult,max_mult,xTitle,"n_{4}") );
    h_Q1_vsMult.push_back(       createProfile(createName(bn,pfn,"Q1"),      nBins_mult,min_mult,max_mult,xTitle,"Q_{1}") );
    h_Q1Q1_vsMult.push_back(     createProfile(createName(bn,pfn,"Q1Q1"),    nBins_mult,min_mult,max_mult,xTitle,"Q_{1}^{2}") );
    h_Q1Q1Q1_vsMult.push_back(   createProfile(createName(bn,pfn,"Q1Q1Q1"),  nBins_mult,min_mult,max_mult,xTitle,"Q_{1}^{3}") );
    h_Q1Q1Q1Q1_vsMult.push_back( createProfile(createName(bn,pfn,"Q1Q1Q1Q1"),nBins_mult,min_mult,max_mult,xTitle,"Q_{1}^{4}") );
    h_Q2_vsMult.push_back(       createProfile(createName(bn,pfn,"Q2"),      nBins_mult,min_mult,max_mult,xTitle,"Q_{2}^{2}") );
    h_Q2Q1_vsMult.push_back(     createProfile(createName(bn,pfn,"Q2Q1"),    nBins_mult,min_mult,max_mult,xTitle,"Q_{2}Q_{2}") );
    h_Q2Q1Q1_vsMult.push_back(   createProfile(createName(bn,pfn,"Q2Q1Q1"),  nBins_mult,min_mult,max_mult,xTitle,"Q_{2}Q_{1}^2") );
    h_Q2Q2_vsMult.push_back(     createProfile(createName(bn,pfn,"Q2Q2"),    nBins_mult,min_mult,max_mult,xTitle,"Q_{2}^2") );
    h_Q3_vsMult.push_back(       createProfile(createName(bn,pfn,"Q3"),      nBins_mult,min_mult,max_mult,xTitle,"Q_{3}") );
    h_Q3Q1_vsMult.push_back(     createProfile(createName(bn,pfn,"Q3Q1"),    nBins_mult,min_mult,max_mult,xTitle,"Q_{3}Q_{1}") );
    h_Q4_vsMult.push_back(       createProfile(createName(bn,pfn,"Q4"),      nBins_mult,min_mult,max_mult,xTitle,"Q_{4}") );


    }
  if (reportEnd(__FUNCTION__))  {/* no ops */}
}

//________________________________________________________________________
void PtPtHistos::importHistograms(TFile & inputFile)
{
  if (reportStart(__FUNCTION__)) {/* no ops */}
  const String & bn  = getName();
  const String & ptn = getParentName();
  const String & ppn = getParentPathName();
  multName = configuration.getValueInt(ppn,"MultName");
  if (reportInfo(__FUNCTION__))
    {
    printCR();
    printValue("Parent Task Name",ptn);
    printValue("Parent Path Name",ppn);
    printValue("Histo Base Name.",bn);
    printValue("MultName",multName);
    }

  for (auto & particleFilter1 : particleFilters)
    {
    String pfn = particleFilter1->getName();
    h_f1_vsMult.push_back( loadProfile(inputFile,createName(bn,pfn,"n1") ));
    h_f2_vsMult.push_back( loadProfile(inputFile,createName(bn,pfn,"n2") ));
    h_f3_vsMult.push_back( loadProfile(inputFile,createName(bn,pfn,"n3") ));
    h_f4_vsMult.push_back( loadProfile(inputFile,createName(bn,pfn,"n4") ));
    h_Q1_vsMult.push_back(       loadProfile(inputFile,createName(bn,pfn,"Q1") ));
    h_Q1Q1_vsMult.push_back(     loadProfile(inputFile,createName(bn,pfn,"Q1Q1") ));
    h_Q1Q1Q1_vsMult.push_back(   loadProfile(inputFile,createName(bn,pfn,"Q1Q1Q1") ));
    h_Q1Q1Q1Q1_vsMult.push_back( loadProfile(inputFile,createName(bn,pfn,"Q1Q1Q1Q1") ));
    h_Q2_vsMult.push_back(       loadProfile(inputFile,createName(bn,pfn,"Q2") ));
    h_Q2Q1_vsMult.push_back(     loadProfile(inputFile,createName(bn,pfn,"Q2Q1") ));
    h_Q2Q1Q1_vsMult.push_back(   loadProfile(inputFile,createName(bn,pfn,"Q2Q1Q1") ));
    h_Q2Q2_vsMult.push_back(     loadProfile(inputFile,createName(bn,pfn,"Q2Q2") ));
    h_Q3_vsMult.push_back(       loadProfile(inputFile,createName(bn,pfn,"Q3") ));
    h_Q3Q1_vsMult.push_back(     loadProfile(inputFile,createName(bn,pfn,"Q3Q1") ));
    h_Q4_vsMult.push_back(       loadProfile(inputFile,createName(bn,pfn,"Q4") ));
   }
  if (reportEnd(__FUNCTION__))  { /* no ops */ }
}

/// Implement Eqs. D4 - D7
///
void PtPtHistos::fill(double mult,
                      vector<double> & counts,
                      vector<double> & q1Sum,
                      vector<double> & q2Sum,
                      vector<double> & q3Sum,
                      vector<double> & q4Sum)
{
  unsigned int nFilters = counts.size();
  ////unsigned int index = 0;
  for (unsigned int i1=0; i1<nFilters; i1++)
    {
    double n1 = counts[i1];
    double n2 = n1*(n1-1);
    double n3 = n2*(n1-2);
    double n4 = n3*(n1-3);
    h_f1_vsMult[i1]->Fill(mult,n1);
    h_f2_vsMult[i1]->Fill(mult,n2);
    h_f3_vsMult[i1]->Fill(mult,n3);
    h_f4_vsMult[i1]->Fill(mult,n4);

    double Q1 = q1Sum[i1];
    double Q2 = q2Sum[i1];
    double Q3 = q3Sum[i1];
    double Q4 = q4Sum[i1];
    double Q1Q1     = Q1*Q1;
    double Q1Q1Q1   = Q1Q1*Q1;
    double Q1Q1Q1Q1 = Q1Q1Q1*Q1;
    double Q2Q1   = Q2*Q1;
    double Q2Q2   = Q2*Q2;
    double Q2Q1Q1 = Q2Q1*Q1;
    double Q3Q1   = Q3*Q1;

    h_Q1_vsMult[i1]->Fill(mult,Q1);
    h_Q1Q1_vsMult[i1]->Fill(mult,Q1Q1);
    h_Q1Q1Q1_vsMult[i1]->Fill(mult,Q1Q1Q1);
    h_Q1Q1Q1Q1_vsMult[i1]->Fill(mult,Q1Q1Q1Q1);
    h_Q2_vsMult[i1]->Fill(mult,Q2);
    h_Q2Q1_vsMult[i1]->Fill(mult,Q2Q1);
    h_Q2Q1Q1_vsMult[i1]->Fill(mult,Q2Q1Q1);
    h_Q2Q2_vsMult[i1]->Fill(mult,Q2Q2);
    h_Q3_vsMult[i1]->Fill(mult,Q3);
    h_Q3Q1_vsMult[i1]->Fill(mult,Q3Q1);
    h_Q4_vsMult[i1]->Fill(mult,Q4);
    }
}

void PtPtHistos::setParticleFilters(const vector<ParticleFilter*> & _particleFilters)
{
  particleFilters = _particleFilters;
}

vector<ParticleFilter*> & PtPtHistos::getParticleFilters()
{
  return particleFilters;
}

const vector<ParticleFilter*> & PtPtHistos::getParticleFilters() const
{
  return particleFilters;
}


} // namespace CAP

