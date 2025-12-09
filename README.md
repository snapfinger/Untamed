## Untamed: Unconstrained Tensor Decomposition and Graph Node Embedding for Cortical Parcellation

🧠 Introduction
---
We present **Untamed**, a novel framework that integrates unconstrained tensor decomposition using NASCAR to identify functional networks, with state-of-the-art graph node embedding to generate cortical parcellations. Our method produces near-homogeneous, spatially coherent regions aligned with large-scale functional networks, while avoiding strong assumptions like statistical independence required in ICA. Across multiple datasets, Untamed consistently demonstrates improved or comparable performance in functional connectivity homogeneity and task contrast alignment compared to existing atlases. The pipeline is fully automated, allowing for rapid adaptation to new datasets and the generation of custom parcellations. The atlases derived from the Genomics Superstruct Project (GSP) dataset, along with the code for generating customizable parcel numbers, are publicly available at:

👉 https://untamed-atlas.github.io

🧩 Pipeline Overview
---
![Untamed Framework](./figs/untamed_framework.jpg)

The Untamed pipeline consists of three main stages:
1. **Graph Construction**: Build adjacency matrices from NASCAR spatial maps with spatial constraints
2. **Graph Embedding**: Generate node embeddings using NetMF (Network Matrix Factorization)
3. **Clustering**: Apply k-means clustering to produce the final parcellation

---

## 📋 Table of Contents

- [Installation & Requirements](#-installation--requirements)
- [Quick Start Guide](#-quick-start-guide)
- [Step-by-Step Tutorial](#-step-by-step-tutorial)
- [Using with Custom Data](#-using-with-custom-data)
- [Parameter Guide](#%EF%B8%8F-parameter-guide)
- [Output Description](#-output-description)
- [Troubleshooting](#-troubleshooting)
- [Citation](#-citation)

---

## 🔧 Installation & Requirements

### System Requirements
- **MATLAB** (R2020b or later recommended)
- Sufficient memory (100GB+ RAM recommended)

### Installation Steps

1. **Clone or download this repository**
   ```bash
   git clone https://github.com/snapfinger/Untamed.git
   cd Untamed
   ```

2. **Verify MATLAB installation**
   - Open MATLAB and ensure it's properly configured
   - All required functions are included in the repository

3. **Optional Dependencies** (for generating NASCAR networks from scratch)
   - **TFOCS toolbox** (http://cvxr.com/tfocs/): Required only if you need to generate NASCAR networks with constraints
   - **N-way Toolbox** (http://www.models.life.ku.dk/nwaytoolbox/download): Required only for CORCONDIA metric calculation

   *Note: These are only needed if you're generating NASCAR networks from raw fMRI data. For using pre-computed networks (as in the quick start), they are not required.*

---

## 🚀 Quick Start Guide

### Using Pre-computed NASCAR Networks (Recommended for First-Time Users)

1. **Navigate to the repository directory in MATLAB**
   ```matlab
   cd /path/to/Untamed
   ```

2. **Run the main script**
   ```matlab
   run_untamed.m
   ```

3. **Find your results**
   - Outputs will be saved in `./untamed_output/`
   - The final parcellation is in `parcel.mat`

That's it! The script will automatically:
- Load the pre-computed NASCAR networks from `artifact/nascar_GSP_rst.mat`
- Construct graph adjacency matrices
- Generate embeddings
- Create a parcellation with 50 parcels per hemisphere (default)

---

## 📖 Step-by-Step Tutorial

This section provides a detailed walkthrough of the Untamed pipeline, explaining each step and how to customize it for your needs.

### Step 1: Understanding the Input Data

The pipeline requires **NASCAR spatial maps** as input. These are functional network spatial patterns derived from tensor decomposition of multi-subject fMRI data.

**Pre-computed input (included):**
- `artifact/nascar_GSP_rst.mat`: Contains NASCAR networks from the GSP dataset
  - Structure: `rst_50{2, 1}` contains the spatial maps
  - Format: `[num_vertices × num_networks]` matrix
  - Surface space: fsaverage6 (37,476 left + 37,471 right = 74,947 vertices)

**What if you need to generate NASCAR networks from your own data?**
- See the [Using with Custom Data](#using-with-custom-data) section below
- Refer to the NASCAR documentation in `helper/NASCAR_v2/README`
- Example code is available in `helper/NASCAR_v2/demo.m`

### Step 2: Graph Construction

The first stage builds a graph representation of cortical connectivity.

**What happens:**
1. **Correlation computation**: Spatial correlation matrix is computed from NASCAR maps
   ```matlab
   spatial = rst_50{2, 1};
   nascar_corr = corr(spatial');
   ```

2. **Hemisphere separation**: Left and right hemispheres are processed separately
   ```matlab
   nascar_corr_L = nascar_corr(1:37476, 1:37476);
   nascar_corr_R = nascar_corr(37477:end, 37477:end);
   ```

3. **Gaussian kernel transformation**: Converts correlations to similarity weights
   ```matlab
   sigma = 0.5;
   nascar_corr_L = exp(nascar_corr_L ./ (2 * sigma^2));
   ```

4. **Spatial constraint addition**: Ensures spatial coherence by restricting connections to neighboring vertices
   ```matlab
   nascar_corrL_nb = addSpatialConstraintOneHem(nascar_corr_L, 'L', nb_L);
   ```

**Output:** Adjacency matrices saved to `./untamed_output/adj.mat`

### Step 3: Graph Embedding

The second stage generates low-dimensional embeddings for each cortical vertex.

**What happens:**
1. **NetMF embedding**: Applies Network Matrix Factorization to learn vertex representations
   ```matlab
   embed_L = netmf_embed_wrapper(nascar_corrL_nb, negative_sampling, alpha, T, embed_dim);
   ```

2. **Parameters:**
   - `T = 7`: Window length for random walk approximation
   - `embed_dim = 128`: Dimensionality of output embeddings
   - `alpha = 0.5`: Power scaling factor
   - `negative_sampling = 1`: Negative sampling flag

**Output:** Embeddings saved to `./untamed_output/embed.mat`

### Step 4: Clustering

The final stage clusters embeddings to produce the parcellation.

**What happens:**
1. **Multiple k-means runs**: Performs 500 independent k-means clustering runs with random initializations
   ```matlab
   [parcels, sumd_exps] = kmeans_cluster(embed_L, embed_R, ...
                                         cluster_num_left, cluster_num_right, ...
                                         km_max_iter, exp_num);
   ```

2. **Best parcellation selection**: Chooses the run with the smallest within-cluster sum of squares
   ```matlab
   [parcel, ~] = kmeans_best_cluster(parcels, sumd_exps);
   ```

**Output:** Final parcellation saved to `./untamed_output/parcel.mat`

---

## 🔬 Using with Custom Data

### Option 1: Using Your Own Brain Networks

If you have pre-computed brain networks:

1. **Prepare your brain network file:**
   - Format: MATLAB `.mat` file
   - Variable: `rst_50` (or adjust the variable name in `run_untamed.m`)
   - Structure: `rst_50{2, 1}` should contain `[num_vertices × num_networks]` spatial maps
   - Surface space: Must be in fsaverage6 space (or adjust vertex numbers accordingly)

2. **Update the script:**
   ```matlab
   % In run_untamed.m, line 27, change:
   load('artifact/nascar_GSP_rst.mat');
   % To:
   load('path/to/your/brain_network_file.mat');
   ```

3. **Adjust vertex numbers** (if not using fsaverage6):
   ```matlab
   % In run_untamed.m, lines 33-34, update:
   LEFT_VERT_NUM = your_left_vertex_count;
   RIGHT_VERT_NUM = your_right_vertex_count;
   ```

4. **Update spatial constraint files** (if needed):
   - The script uses `files/fsavg6_nb1_index.mat` for spatial constraints
   - For different surface spaces, you'll need to generate corresponding neighborhood files

### Option 2: Generating NASCAR Networks from Raw fMRI Data

If you need to generate NASCAR networks from scratch:

1. **Prepare your fMRI data:**
   - Format: 3D tensor `[vertices × timepoints × subjects]`
   - Preprocessing: Standard fMRI preprocessing (motion correction, registration to fsaverage6, etc.)

2. **Temporal alignment (recommended especially for non-task fMRI data):**
   ```matlab
   % Use GroupBrainSync for temporal alignment
   addpath('helper/GroupBrainSync');
   aligned_data = groupBrainSync(your_fmri_tensor);
   ```

3. **Run NASCAR decomposition:**
   ```matlab
   addpath('helper/NASCAR_v2');
   option = srscpd('opt');
   option.nonnegative = [0 0 1];  % Non-negative constraint on spatial mode
   option.rankOneOptALS.useTFOCS = false;
   option.optAlg.normRegParam = 0.001;
   option.optAlg.optSolver.maxNumIter = 2000;
   
   R = 50;  % Desired number of networks
   result = srscpd(aligned_data, R, option);
   
   % Extract spatial maps
   spatial_maps = result(R).U{1};  % [vertices × R]
   ```

4. **Save in the expected format:**
   ```matlab
   rst_50{2, 1} = spatial_maps;
   save('your_nascar_output.mat', 'rst_50');
   ```

5. **Use the saved file** as described in Option 1 above.

---

## ⚙️ Parameter Guide

### Key Hyperparameters

These parameters control the behavior of the Untamed pipeline and can be adjusted in `run_untamed.m`:

#### Graph Construction Parameters

- **`nb`** (default: `1` for testing purposes only, use 20-50 for real data)
  - Spatial neighborhood constraint level
  - Controls how many neighbors are considered for spatial constraints
  - Higher values allow connections to more distant neighbors

- **`sigma`** (default: `0.5`)
  - Gaussian kernel bandwidth for similarity transformation
  - Controls the decay rate of similarity with correlation
  - Lower values create more sparse graphs
  - Range: 0.1-1.0

#### Embedding Parameters

- **`T`** (default: `7`)
  - NetMF embedding window length
  - Controls the context window for random walk approximation
  - Higher values capture longer-range dependencies
  - Range: 3-10

- **`embed_dim`** (default: `128`)
  - Dimensionality of output embeddings
  - Higher dimensions capture more information but require more memory
  - Range: 128 / 256

- **`alpha`** (default: `0.5`)
  - Power scaling factor for embeddings
  - Controls the influence of singular values
  - Range: 0.0-1.0

- **`negative_sampling`** (default: `1`)
  - Flag for negative sampling in NetMF

#### Clustering Parameters

- **`cluster_num_left`** (default: `50`)
  - Desired number of parcels in the left hemisphere
  - Adjust based on your resolution needs
  - Range: 10-200+

- **`cluster_num_right`** (default: `50`)
  - Desired number of parcels in the right hemisphere
  - Can be different from left hemisphere

- **`km_max_iter`** (default: `20000`)
  - Maximum iterations for each k-means run
  - Increase if convergence issues occur

- **`exp_num`** (default: `500`)
  - Number of independent k-means runs
  - Higher values improve robustness but increase computation time
  - Range: 100-1000

### Parameter Selection Tips

- **For higher resolution parcellations**: Increase `cluster_num_left` and `cluster_num_right`
- **For faster computation**: Reduce `exp_num`, `embed_dim`, or `T`
- **For better spatial coherence**: Increase `nb` or adjust `sigma`
- **For different surface spaces**: Update `LEFT_VERT_NUM` and `RIGHT_VERT_NUM`, and provide corresponding neighborhood files

---

## 📊 Output Description

### Output Files

All outputs are saved in `./untamed_output/` directory:

1. **`adj.mat`**
   - Contains: `nascar_corrL_nb`, `nascar_corrR_nb`
   - Description: Graph adjacency matrices for left and right hemispheres
   - Format: Sparse or full matrices `[num_vertices × num_vertices]`

2. **`embed.mat`**
   - Contains: `embed_L`, `embed_R`
   - Description: Vertex embeddings for left and right hemispheres
   - Format: Cell arrays with embeddings `[num_vertices × embed_dim]`
   - Structure: `embed_L{1, 2}` contains the actual embedding matrix

3. **`parcel.mat`**
   - Contains: `parcel`
   - Description: Final parcellation labels
   - Format: Integer vector `[74947 × 1]` (for fsaverage6)
   - Values: Parcel IDs (1 to `cluster_num_left` for left, `cluster_num_left+1` to `cluster_num_left+cluster_num_right` for right)

### Using the Output

**Load and visualize the parcellation:**
```matlab
load('./untamed_output/parcel.mat');
% parcel is a [74947 × 1] vector with parcel labels
% First 37476 entries: left hemisphere
% Remaining entries: right hemisphere
```

**Export for visualization tools:**
```matlab
% For FreeSurfer/Connectome Workbench
% Save as annotation file or gifti format
% (You may need additional scripts for format conversion)
```

---

## 🔍 Troubleshooting

### Common Issues and Solutions

#### Issue: "Out of memory" error
**Solution:**
- Reduce `embed_dim` (e.g., from 128 to 64)
- Process hemispheres separately if already not doing so
- Close other MATLAB processes
- Use a machine with more RAM

#### Issue: K-means not converging
**Solution:**
- Increase `km_max_iter` (e.g., to 50000)
- Check if embeddings contain NaN or Inf values
- Try different random seeds

#### Issue: Parcellation looks fragmented
**Solution:**
- Increase `nb` (spatial neighborhood constraint)
- Adjust `sigma` to create stronger spatial connections
- Check that spatial constraint files are correct for your surface space

#### Issue: Wrong number of vertices
**Solution:**
- Verify `LEFT_VERT_NUM` and `RIGHT_VERT_NUM` match your data
- Check that your networks are in the correct surface space
- Ensure hemisphere separation indices are correct

#### Issue: Can't find required files
**Solution:**
- Ensure you're running from the repository root directory
- Check that `files/fsavg6_nb1_index.mat` exists
- Verify `artifact/nascar_GSP_rst.mat` is present (or your custom file)

#### Issue: Slow computation
**Solution:**
- Reduce `exp_num` (fewer k-means runs)
- Reduce `embed_dim` or `T`
- Use a machine with more CPU cores (some functions support parallel processing)
- Check if you can use sparse matrices for adjacency matrices

### Getting Help

- Check the code comments in `run_untamed.m` for detailed explanations
- Review the NASCAR documentation in `helper/NASCAR_v2/README`
- Ensure all required files are in the correct locations
- Verify MATLAB version compatibility (R2016b+)

---

## 📁 Repository Contents

- **`run_untamed.m`**: Main script to generate parcellations
- **`artifact/nascar_GSP_rst.mat`**: Pre-computed NASCAR networks from GSP dataset
- **`helper/`**: Supporting functions and toolboxes
  - `NASCAR_v2/`: NASCAR tensor decomposition implementation
  - `GroupBrainSync/`: Temporal alignment for multi-subject data
  - NetMF embedding functions
  - Clustering utilities
- **`files/`**: Required data files (spatial constraints, vertex mappings)
- **`vis_spatial_corr.ipynb`**: Jupyter notebook for visualizing spatial correlations of NASCAR networks and ICA networks
- **`figs/untamed_framework.jpg`**: Pipeline diagram

---

## 📚 Citation

If you use Untamed in your research, please cite:

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

### Additional Citations

If you use the NASCAR component, please also cite:

```bibtex
@article{li2023identification,
  title={Identification of overlapping and interacting networks reveals intrinsic spatiotemporal organization of the human brain},
  author={Li, Jian and Liu, Yijun and Wisnowski, Jessica L and Leahy, Richard M},
  journal={NeuroImage},
  volume={270},
  pages={119944},
  year={2023},
  publisher={Elsevier}
}
```

