## Fun4All TrackSeeding Macro

This repository contains a **Fun4All macro** that runs **TrackSeeding over the cluster DST** in sPHENIX, with optional **KFParticle reconstruction** and **FastHadronReco** support.

---
```cpp
## For K⁰ₛ → π⁺π⁻

fastreco->setPhotonConv(false);
fastreco->setDecayMass1(0.13957);  // pion mass [GeV/c^2]
fastreco->setDecayMass2(0.13957);  // pion mass [GeV/c^2]

## For ϕ → K⁺K⁻

fastreco->setPhotonConv(false);
fastreco->setDecayMass1(0.49367);  // kaon mass [GeV/c^2]
fastreco->setDecayMass2(0.49367);  // kaon mass [GeV/c^2]

## For γ → e⁺e⁻
fastreco->setPhotonConv(true);
fastreco->setDecayMass1(0.000511); // electron mass [GeV/c^2]
fastreco->setDecayMass2(0.000511); // electron mass [GeV/c^2]

