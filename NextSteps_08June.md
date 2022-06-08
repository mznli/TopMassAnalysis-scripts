#### Next Steps (June 08, 2022) 
##### Main Changes
* There are some bugs fixed in the updated code, so please use it.
   * With this script we expect that the Muons from Semi-Lep B-Decays are selected truely. 
      * Check the nSemiMuBDecays histogram now. Has the overflow (events with 3 Bs decaying semi-muonically) gone?
* Also the usage of TLorentzVector object has been shown, where the W boson and Top quark masses are reconstructed at generator level.
   * Check how the two massive objects are constructed.
##### Starting to look at the reconstructed Muons
* The reconstructed Muon collection is called. 
   * Take note how the Muon pt distribution is filled.
      * It's done in two ways: filling a tree branch or using a histogram
   * Compare the pt distributions for soft muons (from semi-leptonic B decays) and hard muons (coming from hard scattering process)
   * Look at the Muon multiplicity distribution. How do you explain it? 
   * Please make plots for other properties of Muons, for example prepare the following histograms:
      * eta, phi, isolation variables, identification criteria and etc.
      * Hint: check the variables from [nanoAOD-documentation page](https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html#Muon)
##### Drop Output Tree Branches (to reduce file size)
* A `keep_and_drop.txt` file is added. Put it in the `nanoAOD-tools/scripts` directory.
###### , and then run:
    python scripts/nano_postproc.py outDir INPUT-ROOT-FILE -I PhysicsTools.NanoAODTools.postprocessing.modules.SemiLepBDecayTagger SemiLepBDecayTagger -N 100 --bo scripts/keep_and_drop.txt
   * This way, the useless branches are dropped. The output file capacity reduces notably.   
