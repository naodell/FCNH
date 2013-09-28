#include "../interface/TCPrimaryVtx.h"
#include "../interface/TCPrimaryVtxLinkDef.h"

TCPrimaryVtx::TCPrimaryVtx() {
}

TCPrimaryVtx::~TCPrimaryVtx() {
}

float TCPrimaryVtx::NDof() const {
   return _nDof;
}

float TCPrimaryVtx::Chi2() const {
   return _chi2;
}

bool TCPrimaryVtx::IsFake() const {
   return _isFake;
}

int TCPrimaryVtx::Ntracks() const {
   return _nTracks;
}

float TCPrimaryVtx::SumPt2Trks() const {
  return _sumPt2Trks;
}


void TCPrimaryVtx::SetNDof(float n) {
   _nDof = n;
}

void TCPrimaryVtx::SetChi2(float chi2) {
   _chi2 = chi2;
}

void TCPrimaryVtx::SetIsFake(bool isF) {
   _isFake = isF;
}

void TCPrimaryVtx::SetNtracks(int nTrk) {
   _nTracks = nTrk;
}

void TCPrimaryVtx::SetSumPt2Trks(float sumPt2) {
  _sumPt2Trks = sumPt2;
}
