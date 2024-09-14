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
#include "ParticleDecayerTask.hpp"
#include "PrintHelpers.hpp"

ClassImp(CAP::ParticleDecayerTask);

namespace CAP
{

ParticleDecayerTask::ParticleDecayerTask()
:
EventTask(),
decayer()
{
  appendClassName("ParticleDecayerTask");
  setName("ParticleDecayerTask");
  setTitle("ParticleDecayerTask");
  setVersion("1.0");
 }

void ParticleDecayerTask::setDefaultConfiguration()
{
  EventTask::setDefaultConfiguration();
  addProperty("Stream0Index", 0);
}

// ====================================================================
// Decay unstable particles in the event
// Decayed particles are retained in the event but their live
// flag is set to false.
// Children particles are added at the tail end of the event
// with a live flag  set to true.
// ====================================================================
void ParticleDecayerTask::execute()
{
  Event & event = *Manager<Event>::getObjectAt(0);
  int nParticles = event.getNParticles();
  if (nParticles< 1)
    {
    if (reportWarning(__FUNCTION__)) cout << "Skipping empty event." << endl;
    return;
    }

  unsigned int iParticle = 0;
  bool done = false;
  while (!done)
    {
    Particle & parent = * event.getParticleAt(iParticle);
    if (parent.isLive() &&  !parent.isStable() )
      {
      ParticleType   &    parentType      = parent.getType();
      LorentzVector &    parentMomentum  = parent.getMomentum();
      LorentzVector &    parentPosition  = parent.getPosition();
      ParticleDecayMode & decayMode       = parentType.generateDecayMode();
      int nChildren = decayMode.getNChildren();
      switch (nChildren)
        {
          case 1:
          if (reportInfo(__FUNCTION__)) cout << "case 1  parentType==" << parent.getName() << endl;
          break;
          case 2:
          {
          Particle * child1 = particleFactory->getNextObject();
          Particle * child2 = particleFactory->getNextObject();
          ParticleType   & childType1 = decayMode.getChildType(0); child1->setType(&childType1); child1->setLive(true);
          ParticleType   & childType2 = decayMode.getChildType(1); child2->setType(&childType2); child2->setLive(true);
          LorentzVector & p1 = child1->getMomentum();
          LorentzVector & r1 = child1->getPosition();
          LorentzVector & p2 = child2->getMomentum();
          LorentzVector & r2 = child2->getPosition();
          decayer.decay2(parentType,
                         parentMomentum,
                         parentPosition,
                         childType1,p1,r1,
                         childType2,p2,r2);
          event.addParticle(child1); child1->setLive(true);
          event.addParticle(child2); child2->setLive(true);
          parent.setDecayed(true);
          nParticles += 2;
          }
          break;

          case 3:
          {
          Particle * child1 = particleFactory->getNextObject();
          Particle * child2 = particleFactory->getNextObject();
          Particle * child3 = particleFactory->getNextObject();
          ParticleType   & childType1 = decayMode.getChildType(0); child1->setType(&childType1); child1->setLive(true);
          ParticleType   & childType2 = decayMode.getChildType(1); child2->setType(&childType2); child2->setLive(true);
          ParticleType   & childType3 = decayMode.getChildType(2); child3->setType(&childType3); child3->setLive(true);
          LorentzVector & p1 = child1->getMomentum();
          LorentzVector & r1 = child1->getPosition();
          LorentzVector & p2 = child2->getMomentum();
          LorentzVector & r2 = child2->getPosition();
          LorentzVector & p3 = child3->getMomentum();
          LorentzVector & r3 = child3->getPosition();
          decayer.decay3(parentType,
                         parentMomentum,
                         parentPosition,
                         childType1,p1,r1,
                         childType2,p2,r2,
                         childType3,p3,r3);
          event.addParticle(child1); child1->setLive(true);
          event.addParticle(child2); child2->setLive(true);
          event.addParticle(child3); child3->setLive(true);
          parent.setDecayed(true);
          nParticles += 3;
          }
          break;
          case 4:
//          {
//  //        Particle * child1 = particleFactory->getNextObject();
//  //        Particle * child2 = particleFactory->getNextObject();
//  //        Particle * child3 = particleFactory->getNextObject();
//  //        Particle * child4 = particleFactory->getNextObject();
//  //        ParticleType   & childType1 = decayMode.getChildType(0); child1->setType(&childType1); child1->setLive(true);
//  //        ParticleType   & childType2 = decayMode.getChildType(1); child2->setType(&childType2); child2->setLive(true);
//  //        ParticleType   & childType3 = decayMode.getChildType(2); child3->setType(&childType3); child3->setLive(true);
//  //        ParticleType   & childType4 = decayMode.getChildType(3); child4->setType(&childType4); child4->setLive(true);
//  //        LorentzVector & p1 = child1->getMomentum();
//  //        LorentzVector & r1 = child1->getPosition();
//  //        LorentzVector & p2 = child2->getMomentum();
//  //        LorentzVector & r2 = child2->getPosition();
//  //        LorentzVector & p3 = child3->getMomentum();
//  //        LorentzVector & r3 = child3->getPosition();
//  //        LorentzVector & p4 = child4->getMomentum();
//  //        LorentzVector & r4 = child4->getPosition();
//  //        decayer.decay4(parentType,
//  //                       parentMomentum,
//  //                       parentPosition,
//  //                       childType1,p1,r1,
//  //                       childType2,p2,r2,
//  //                       childType3,p3,r3,
//  //                       childType4,p4,r4);
//  //        event.add(child1); child1->setLive(true);
//  //        event.add(child2); child2->setLive(true);
//  //        event.add(child3); child3->setLive(true);
//  //        event.add(child4); child4->setLive(true);
//  //        parent.setDecayed(true);
//  //        nParticles += 4;
//          }
          break;
          case 5:
          break;
        }
      }
    iParticle++;
    done = (iParticle >= event.getNParticles());
    }
  nParticles = event.getNParticles();
}





} // namespace CAP