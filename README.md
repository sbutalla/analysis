# analysis
## Stephen D. Butalla

A repository for various scripts for automating high energy physics analysis tasks.

`sensitivity_scan_v5.py`: Automates the `MadGraph5` (MG5) Monte Carlo event simulation software and `MadAnalysis5` (MA5) event analysis software for the "Dark Sector" 
Beyond Standard Model (BSM) theory. The current model implemented (more to follow) is 
pp-> ZD -> fd1 fd1~ -> fd2 fd2~ mu+ mu- mu+ mu-. The aggregate output of this script
is used for performing sensitivity calculations (production cross-section, branching ratio, and geometric/kinematic acceptances) to determine if the siganture of this
model will be possible to detect in data obtained from the Compact Muon Solenoid (CMS) experiment during Runs 1 and 2, as well as Run 3 and beyond. For a given (or range of) ZD, fd1, and fd2 mass(es), the script uses the Universal FeynRules Output (UFO) files produced from the Mathematic package `FeynRules`. A new model is built in MadGraph from these UFO files, initialized with the mass or masses input at the execution of the script. At this time, the option to calculate the decay widths (and branching ratios) for the decay is performed. Production cross-section of the process pp -> ZD is calculated as a separate process, as the production cross section from the interaction is incorrect. The outut of this new model is then used to generate 10000 events, creating the Les Houches Event (LHE) file for analysis by MA5. The geometric and kinematic acceptance calculation is performed.
