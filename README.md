## Untamed: Unconstrained Tensor Decomposition and Graph Node Embedding for Cortical Parcellation

🧠 Introduction
---
We present **Untamed**, a novel framework that integrates unconstrained tensor decomposition using NASCAR to identify functional networks, with state-of-the-art graph node embedding to generate cortical parcellations. Our method produces near-homogeneous, spatially coherent regions aligned with large-scale functional networks, while avoiding strong assumptions like statistical independence required in ICA. Across multiple datasets, Untamed consistently demonstrates improved or comparable performance in functional connectivity homogeneity and task contrast alignment compared to existing atlases. The pipeline is fully automated, allowing for rapid adaptation to new datasets and the generation of custom parcellations. The atlases derived from the Genomics Superstruct Project (GSP) dataset, along with the code for generating customizable parcel numbers, are publicly available at:

👉 https://untamed-atlas.github.io

🧩 Pipeline
---
![Untamed Framework](./figs/untamed_framework.jpg)


📁 Content
----
- [`artifact/nascar_GSP_rst.mat`](artifact/nascar_GSP_rst.mat): Contains NASCAR networks derived from the [GSP dataset](https://www.nature.com/articles/sdata201531).

- [`vis_spatial_corr.ipynb`](vis_spatial_corr.ipynb): Jupyter notebook for visualizing spatial correlations between brain networks. It compares ICA components from the [HCP S1200 release](https://www.humanconnectome.org/study/hcp-young-adult/document/1200-subjects-data-release) (available on the [HCP website](https://db.humanconnectome.org/app/template/Login.vm)) with the NASCAR networks in the file above.



🚀 Generate your own *Untamed* atlas
---
1. Install MATLAB (if not already installed).

2. Replace `nascar_GSP_rst.mat` with a network file generated from your own dataset.

3. Run `run_untamed.m`` to produce the parcellation.

📚 Citation
---
```bibtex
@article{untamed2025,
    title={Untamed: Unconstrained Tensor Decomposition and Graph Node Embedding for Cortical Parcellation},
    author={Liu, Yijun and Li, Jian and Wisnowski, Jessica L and Leahy, Richard M},
    journal={bioRxiv},
    pages={2024.01.05.574423},
    year={2025},
    publisher={Cold Spring Harbor Laboratory}
}
```

