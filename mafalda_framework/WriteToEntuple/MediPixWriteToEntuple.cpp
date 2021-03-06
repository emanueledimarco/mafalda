/**
 * Author: John Idarraga <idarraga@cern.ch> , 2008
 * Medipix Group, Universite de Montreal
 *
 */

#ifndef MediPixWriteToEntuple_cxx
#define MediPixWriteToEntuple_cxx

#include <vector>
#include <iostream>
#include "MediPixWriteToEntuple.h"


WriteToNtuple::WriteToNtuple(){

  outputROOTFile = new TFile("MPXResults.root","RECREATE");

}

WriteToNtuple::WriteToNtuple(TString fileName){

  outputROOTFile = new TFile(fileName,"RECREATE");

}

void WriteToNtuple::includeTree(TString treeName){

  outputTrees[treeName] = new TTree(treeName, treeName);

}

void WriteToNtuple::writeTree(TString treeName){

  outputTrees[treeName]->Write();

}

void WriteToNtuple::closeNtuple(){

  outputROOTFile->Close();

}

TTree * WriteToNtuple::getTree(TString algoName){

  std::map<TString, TTree *>::iterator treesItr;
  TTree * nothing = 0;
  
  return outputTrees[algoName];

  for ( treesItr = outputTrees.begin() ; treesItr != outputTrees.end(); treesItr++ )
    {
      if((*treesItr).first == algoName)
	{
	  return (*treesItr).second;
	}
    }

  return nothing;
}

/*
Int_t WriteToNtuple::fillVars(std::vector<Int_t> overlappedV_i, std::vector<Int_t> overlappedSumV_i, std::vector<Int_t> nonOverlappedV_i){

  std::vector<Int_t>::iterator it = nonOverlappedV_i.begin();
  Int_t cntr1 = 0;

  if((Int_t)overlappedV_i.size() * (Int_t)nonOverlappedV_i.size() == 0)
    return 0;

  //std::cout << "non Over : " << std::endl;
  for( ; it < nonOverlappedV_i.end(); it++ )
    {
      if(cntr1 < MAX_TO_NTUPLE)
	muons.nonOverlapped[cntr1++] = *it;
    }
  muons.nNonOver = --cntr1;

  // overlapped
  it = overlappedV_i.begin();
  cntr1 = 0;
  for( ; it < overlappedV_i.end(); it++ )
    {
      //std::cout << cntr2++ << ": " << *it << "  , ";
      if(cntr1 < MAX_TO_NTUPLE)
	muons.overlapped[cntr1++] = *it;
    }
  muons.nOver = --cntr1;

  // overlapped sum
  it = overlappedSumV_i.begin();
  cntr1 = 0;
  for( ; it < overlappedSumV_i.end(); it++ )
    {
      //std::cout << cntr2++ << ": " << *it << "  , ";
      if(cntr1 < MAX_TO_NTUPLE)
	{
	  muons.overlappedSum[cntr1++] = *it;
	}
    }
  muons.nOverSum = --cntr1;

  outputTree->Fill();
  //cleanUpArrays();

  return 1;
}





void WriteToNtuple::cleanUpArrays(){

  for(Int_t itr = 0 ; itr < MAX_TO_NTUPLE ; itr++)
    {
      //overlapped[itr] = -1;
      //nonOverlapped[itr] = -1;
    }

}

*/

#endif
