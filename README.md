### Untamed
Code repository to generate Untamed atlas.



#### Abstract
Cortical parcellation has long been a cornerstone in the field of neuroscience, enabling the cerebral cortex to be partitioned into distinct, non-overlapping regions that facilitate the interpretation and comparison of complex neuroscientific data. In recent years, these parcellations have frequently been based on the use of resting-state fMRI (rsfMRI) data. In parallel, methods such as independent components analysis have long been used to identify large-scale functional networks with significant spatial overlap between networks. Despite the fact that both forms of decomposition make use of the same spontaneous brain activity measured with rsfMRI, a gap persists in establishing a clear relationship between disjoint cortical parcellations and brain-wide networks. To address this, we introduce a novel parcellation framework that integrates NASCAR, a three-dimensional tensor decomposition method that identifies a series of functional brain networks, with state-of-the-art graph representation learning to produce cortical parcellations that represent near-homogeneous functional regions that are consistent with these brain networks.  Further, through the use of the tensor decomposition, we avoid the limitations of traditional approaches that assume statistical independence or orthogonality in defining the underlying networks. Our findings demonstrate that these parcellations are comparable or superior to established atlases in terms of homogeneity of the functional connectivity across parcels, task contrast alignment, and architectonic map alignment. Our methodological pipeline is highly automated, allowing for rapid adaptation to new datasets and the generation of custom parcellations in just minutes, a significant advancement over methods that require extensive manual input. We describe this integrated approach, which we refer to as Untamed, as a tool for use in the fields of cognitive and clinical neuroscientific research. Parcellations created from the Human Connectome Project dataset using Untamed, along with the code to generate atlases with custom parcel numbers, are publicly available at https://untamed-atlas.github.io.


#### Citation
```
@article{liu2024graph_bioarxiv,
    title={Graph Learning for Cortical Parcellation from Tensor Decompositions of Resting-State fMRI},
    author={Liu, Yijun and Li, Jian and Wisnowski, Jessica L and Leahy, Richard M},
    journal={bioRxiv},
    pages={2024--01},
    year={2024},
    publisher={Cold Spring Harbor Laboratory}
}
```

