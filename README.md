# TopMassAnalysis-scripts
Simple scripts to prepare for a Top quark mass measurement project
### First you need to prepare the environment and get the nanoAOD tools from the following link:
    https://github.com/cms-nanoAOD/nanoAOD-tools
### Then simply put the script provided in:
    python/postprocessing/modules/
and then run the code:

    python scripts/nano_postproc.py outDir INPUT-ROOT-FILE -I PhysicsTools.NanoAODTools.postprocessing.modules.SemiLepBDecayTagger SemiLepBDecayTagger -N 100
-------------------------------------    
#### Perform some tests
* Check the lepton multiplicity histograms produced from the scripts
* How many events contain a semi-leptonic B decay? Is that reasonable?
* Now increase the statistics and run the code over 10000 events. Check the number of events with more than 1 leptonic B hadron decay (this can be done by looking at nSemiMuBDecays histogram). Check the overflow of the same histogram. How do you explain it?
