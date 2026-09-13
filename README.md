Codes and figures used in the article Detecting epidemic outbreaks with routine blood tests: A comparison of surveillance methods using COVID-19 simulations and real-world seasonal disease data

### Abstract

Timely detection of an epidemic outbreak is essential for implementing interventions and limiting the spread of the pathogen. The digitalization of healthcare has led to the availability of up-to-date, large-scale data on routinely tested blood parameters. This data can be used to detect the onset of an epidemic, which is particularly important when a new, previously unknown pathogen emerges for which no testing systems are available. In this study, we evaluated the effectiveness of various mathematical approaches for detecting an epidemic based on routinely collected blood test results. The evaluation was performed using both simulated and real-world data. First, we developed the BloodCovasim agent-based model, which defines the dependency of C-reactive protein (CRP) on time for each individual based on the course of their disease and simulates blood analysis. Using synthetic data from thousands of BloodCovasim simulations, we determined the characteristics of different methods of epidemic outbreak detections under controlled conditions. Second, we used the same methods but applied them to real-world data on CRP from Russia during a period of increased seasonal infectious diseases (rhinoviruses, adenoviruses, seasonal coronaviruses) associated with the back-to-school effect. Third, we compared the CRP-based detection method with methods based on a set of inflammatory biomarkers obtained from complete blood count (CBC) results. As a result, we found that the applicability of the methods depends on the volume of incoming data. With large amounts of data, they demonstrate similar results in terms of the probability of detecting an outbreak. As the number of tests decreases, detectors based on basic statistical tests lose accuracy faster than those based on moving averages, while the AIC-based detector performs best. The use of composite biomarkers increases the probability of detecting an outbreak compared to using CRP alone. These results can help epidemiologists select mathematical approaches based on the volume of incoming data.

### Repository structure

```
.
├── figures
│   ├── figure 1 (abstract)
│   │   ├── graph
│   │   └── icons
│   ├── figure 2 (bloodcovasim)
│   │   ├── panel A (crp_dynamics)
│   │   │   └── icons
│   │   └── panel B (flowchart)
│   ├── figure 3 (aic)
│   │   └── graph
│   ├── figure 4 (synthetic)
│   ├── figure 5 (real-world)
│   ├── figure 6 (biomarkers)
│   └── supplementary
│       ├── s1 figure (synthetic_figure)
│       ├── s2 figure (real_figure)
│       └── s3 figure (biomarkers)
├── references
└── scripts
```

### Citation
If you find this work helpful and would like to use our model and data, please cite us with
```
@article{
  $add
}

```
### Contact us
Please contact the corresponding author (email: konstantin.a.klochkov (a) gmail.com) if you have any questions.
