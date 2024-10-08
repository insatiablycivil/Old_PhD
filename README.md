# Trip Current Trace Analysis of Medium Voltage Circuitbreakers (Old PhD)

## Scope
My PhD was initially (2015-2019) looking at trip-coil current traces from medium voltage circuitbreakers. 
Goal was to apply "data analytics". After extensive desk study, I manually created labelled dataset using industry data.
This process also allowed me to give insights into common issues in the measuring equipment, and to provide updated guidelines to the industrial partner regarding "prototypical" traces that I synthetically generated and timings for different circuitbreaker models, which they adopted.
My initial work looked to using this labelled dataset for two purposes:
  1. Automated validation of existing automated feature extraction process ([Link to draft of conference paper](https://pureportal.strath.ac.uk/en/publications/automated-feature-validation-of-trip-coil-analysis-in-condition-m))
  2. Automated diagnosis partially using existing automated feature extraction process ([Link to draft of conference paper](https://strathprints.strath.ac.uk/64784/1/Hosseini_etal_ISGTC_2018_Current_based_trip_coil_analysis_of_circuit_breakers.pdf))

The code segments the traces based on the existing automated feature extraction process, then generates new features based on my understanding of what common failure modes are and how they appear symptomatically within the trace. These features then feed into an ensemble of now very basic machine learning approaches (decision trees, svms, etc). At the time, Matlab did not have the functionality to do this more elegantly, perhaps it does now. I no longer have Matlab so cannot actually run this code to check it works / is correct version. In many ways, I find the code quite embarrassing, but at the same time, that perhaps indicates I have grown as a programmer. 

I was working on something much more elaborate (and probably out of my depth), but I took a break once my father died and found it difficult to return to this topic. I instead focussed on the Transformer DGA, and ultimately wrote a new thesis on that topic. It involved a physics-based model using SimScape in Matlab that would attempt to parameter-tune itself to fit a given trace. Its ability to match the trace, as measured via Dynamic Time Warping was to be the core metric for anomaly detection. It was very much WiP. Since I can't run anything, hard to know where it is; I have uploaded three drafts (no idea what is in them sorry!) I have found in case it is useful to anyone.
