## Untamed: Unconstrained Tensor Decomposition and Graph Node Embedding for Cortical Parcellation


Cortical parcellation is fundamental to neuroscience, enabling the division of cerebral cortex into distinct, non-overlapping regions to support interpretation and comparison of complex neuroimaging data. Although extensive literature has investigated cortical parcellation and its connection to functional brain networks, the optimal spatial features for deriving parcellations from resting-state fMRI (rsfMRI) remain unclear. Traditional methods such as Independent Component Analysis (ICA) have been widely used to identify large-scale functional networks, while other approaches define disjoint cortical parcellations. However, bridging these perspectives through effective feature extraction remains an open challenge. To address this, we introduce **Untamed*, a novel framework that integrates unconstrained tensor decomposition using NASCAR to identify functional networks, with state-of-the-art graph node embedding to generate cortical parcellations. Our method produces near-homogeneous, spatially coherent regions aligned with large-scale functional networks, while avoiding strong assumptions like statistical independence required in ICA. Across multiple datasets, Untamed consistently demonstrates improved or comparable performance in functional connectivity homogeneity and task contrast alignment compared to existing atlases. The pipeline is fully automated, allowing for rapid adaptation to new datasets and the generation of custom parcellations. The atlases derived from the Genomics Superstruct Project (GSP) dataset, along with the code for generating customizable parcel numbers, are publicly available at https://untamed-atlas.github.io.


![Untamed Framework](./figs/untamed_framework.jpg)

#### Citation
```
@article{untamed2025,
    title={Untamed: Unconstrained Tensor Decomposition and Graph Node Embedding for Cortical Parcellation},
    author={Liu, Yijun and Li, Jian and Wisnowski, Jessica L and Leahy, Richard M},
    journal={bioRxiv},
    pages={2024.01.05.574423},
    year={2025},
    publisher={Cold Spring Harbor Laboratory}
}
```

