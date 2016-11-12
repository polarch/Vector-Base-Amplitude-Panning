# Vector Base Amplitude Panning Library
#### A compact library implementing the Vector Base Amplitude Panning (VBAP) method and variants for sound spatialization.

---
>
>    Archontis Politis, 2015  
>
>    Department of Signal Processing and Acoustics, Aalto University, Finland  
>
>    archontis.politis@aalto.fi
>
---

This Matlab/Octave library was developed during my doctoral research in the [Communication Acoustics Research Group] (http://spa.aalto.fi/en/research/research_groups/communication_acoustics/), Aalto University, Finland. If you would like to reference the code, you can refer to my dissertation published [here](https://aaltodoc.aalto.fi/handle/123456789/22499):

    Archontis Politis, Microphone array processing for parametric spatial audio techniques, 2016
    Doctoral Dissertation, Department of Signal Processing and Acoustics, Aalto University, Finland
    
## Description

This is a compact Matlab/Octave library implementing vector-base amplitude
panning (VBAP) [ref.1], VBAP-base spreading of a panned source [ref.2&3], 
and Multiple-direction amplitude panning (MDAP) [ref.3]. A function 
implementing the variant vector-base intensity panning [ref.4] is also 
included. Recently, both VBAP and VBIP have been additionaly used for the 
design of robust ambisonic decoding matrices, see [ref.5&6].

The code is written by Archontis Politis, except the core vbap()
function contributed by Ville Pulkki, with small modifications by
Archontis Politis. The following code examples are meant to give a quick
idea how to use the library for common operations in amplitude panning,
such as triangulation of a 3D loudspeaker setup into loudspeaker triplets, 
spreading of a panned source, construction of panning gain
tables, and panning a moving source in a real-time block processing context.

For a demonstration of the various function, see [http://research.spa.aalto.fi/projects/vbap-lib/vbap.html],
or go through the TEST_VBAP_SCRIPTS.m script. For more detailed information, check the help output of the functions 
and the code itself.

The library contains the following main functions:
  
* findLsPairs:    find sorted loudspeaker pairs from loudspeaker directions (for 2D layouts)
* findLSTriplets: find valid loudspeaker triangles from loudspeaker directions (for 3D layouts)
* invertLsMtx:    precompute inversion of matrix of loudspeaker triplets or pairs, for use in VBAP
* getSpreadSrcDirs:   get auxiliary source directions around panning direction, for source spreading and MDAP
* vbap:   Return VBAP panning gains for multiple panning directions, with spread control if needed

Additionaly:

* plotTriangulation:  Plots the loudspeaker triangulated mesh
* getGainTable:   Construct a look-up VBAP gain table of VBAP for a specified regular grid
* vbip:   Similar to VBAP, but implementing its energy-based variant (see [ref.4])
* getPValueResponse:  Returns VBAP frequency-dependent normalization values, for approximate flat perceived response of a panned source in dry playback environments (see [ref.7])

For any questions, comments, corrections, or general feedback, please
contact archontis.politis@aalto.fi

## References

    [1] Pulkki, V. (1997). Virtual Sound Source Positioning Using Vector Base Amplitude Panning. 
    Journal of the Audio Engineering Society, 45(6), 456-466.

    [2] Pulkki, V. (2000). Generic panning tools for MAX/MSP.
    International Computer Music Conference (ICMC), Berlin, Germany

    [3] Pulkki, V. (1999). Uniform Spreading of Amplitude Panned Sources. 
    IEEE Workshop on Applications of Signal Processing to Audio and Acoustics (WASPAA), New Paltz, NY, USA

    [4] Jot, J.-M., Larcher V., Pernaux, J.-M. (1999). A comparative study of 3-D audio encoding and rendering techniques.
    16th International Conference of the AES, Rovaniemi, Finland

    [5] Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. 
    Journal of the Audio Engineering Society, 60(10), 807-820.

    [6] Epain, N., Jin, C.T., Zotter, F. (2014). Ambisonic Decoding With Constant Angular Spread.
    Acta Acustica united with Acustica, 100(May), 928-936.

    [7] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V. (2014). 
    Gain normalization in amplitude panning as a function of frequency and room reverberance. 
    55th International Conference of the AES. Helsinki, Finland.

