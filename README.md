# IntegrateRigor

**A principled framework for gene selection and integration optimization in single-cell data.**

IntegrateRigor brings statistical rigor to two coupled problems in single-cell RNA-seq batch
integration:

1. **Gene selection** — identify *batch-stable* genes (genes whose expression distribution is
   consistent across batches) and separate them from batch-sensitive genes.
2. **Integration optimization** — score the quality of an integrated embedding and automatically
   tune the hyperparameters of an integration method.

Both rely on the same core idea: pick a **reference batch**, fit a mixture model to it, and measure
how much each **query batch** diverges from the reference's structure beyond a simple shift in
mixture proportions. 

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("zhiqianZ/IntegrateRigor")
```

### Dependencies

Core (installed automatically): `Seurat`, `SeuratObject`, `mclust`, `mvtnorm`, `pbmcapply`,
`elbow`, `rlang`, `Matrix`.

Optional, depending on which integration backend you use:

| Backend            | Extra package(s) required                          |
| ------------------ | -------------------------------------------------- |
| CCA / RPCA / Harmony | (none beyond Seurat)                             |
| FastMNN            | `SeuratWrappers`                                    |
| LIGER              | `rliger`                                            |
| scVI               | `reticulate`, `SeuratWrappers`, and a conda env named `scvi` with `scanpy` + `scvi-tools` |

> **Note:** IntegrateRigor expects a Seurat v5 object with a `batch` column (or similarly named)
> in `seurat.obj@meta.data` and raw counts available in the `counts` layer.

---

## Quick start

The end-to-end recommended workflow is:

```r
library(Seurat)
library(IntegrateRigor)

# `obj` is a Seurat v5 object with raw counts and a metadata column "batch".

## 1. Find batch-stable genes (operates on raw counts) --------------------
obj <- BatchStabilityEst(obj, batch = "batch")   # per-gene stability score
obj <- BatchStableGenes(obj)                      # elbow rule -> stable / sensitive sets

stable.genes <- obj@misc$batch.stable.genes

## 2. Preprocess using only batch-stable genes ----------------------------
obj <- Preprocess(obj, batch = "batch", genes = stable.genes, ndims = 30)

## 3. Build a parameter grid and search -----------------------------------
param.df <- make.parameter.df(theta = c(1, 2, 4))      # Harmony's theta, for example

obj <- IntegrateRigor.ParameterS(
  obj,
  parameter.df          = param.df,
  method                = Harmony,
  batch                 = "batch",
  use.batch.stable.genes = TRUE,
  K = 10
)

## 4. Inspect results -----------------------------------------------------
obj@misc$parameter.search                         # scores for every parameter combo
Embeddings(obj, "integrated.bsg.optimal.harmony") # the winning embedding
```

That's the whole loop: **select genes → preprocess → search → use the optimal reduction.** The
sections below explain each step and the alternatives.

---

## Step 1 — Batch-stable gene selection

Genes that behave differently across batches inject technical noise into integration. IntegrateRigor
scores each gene by fitting a **negative-binomial mixture model** (with library size as an offset)
to the reference batch, then measuring how much each query batch deviates from it.

```r
obj <- BatchStabilityEst(
  obj,
  batch     = "batch",   # metadata column holding batch labels
  ref.batch = NULL,       # NULL -> auto-selected via FindReference()
  ngenes    = NULL,       # NULL -> score all sufficiently-expressed genes
  K         = 5,          # number of mixture components
  n.cores   = 10,         # parallel workers (pbmcapply)
  subsample = NULL         # e.g. 0.2 to score on a 20% cell subsample for speed
)
```

Key behaviors:

- If `ref.batch = NULL`, a reference batch is chosen automatically (the batch that covers the most
  clusters) and cached in `obj@misc$reference` for reuse.
- By default, only genes expressed in **> 0.5% of cells** are scored. Pass `ngenes` to instead use
  the top *N* integration features, or `genes` to score a specific set.
- For large datasets (`> 20,000` cells or `> 6` batches), use `subsample` to accelerate.

The per-gene score is written into the RNA assay metadata as `batch.stability.score`. Then apply the
elbow rule to split genes:

```r
obj <- BatchStableGenes(obj, plot = TRUE)   # plot = TRUE shows the elbow cut

obj@misc$batch.stable.genes     # character vector of stable genes
obj@misc$batch.unstable.genes   # batch-sensitive genes
obj@misc$batch.stability.score  # full named score vector
```

---

## Step 2 — Preprocessing

Preprocessing splits the object by batch, normalizes, selects variable features, scales, and runs
PCA — producing the `pca` reduction that the PCA-based integration methods build on.

```r
# For CCA / RPCA / Harmony / FastMNN:
obj <- Preprocess(
  obj,
  batch        = "batch",
  ndims        = 30,
  genes        = stable.genes,  # restrict to batch-stable genes (recommended)
  ngenes       = 2000            # number of variable features to keep
)
```

Backend-specific variants prepare the object differently (no PCA, different normalization):

```r
obj <- Preprocess.scVI(obj, batch = "batch", genes = stable.genes)   # for scVI
obj <- Preprocess.LIGER(obj, batch = "batch", genes = stable.genes)  # for LIGER
```

> **Tip:** To actually integrate *on* batch-stable genes, pass them via `genes =` here. The
> `use.batch.stable.genes` flag in the search functions only controls the naming of the output
> reduction — gene restriction comes from preprocessing.

---

## Step 3 — Integration and scoring (single run)

If you just want to run one integration and score it (rather than a full parameter sweep):

```r
## Run an integration backend
obj <- Integration(
  obj,
  batch        = "batch",
  method       = Harmony,                 # or CCA, RPCA, FastMNN, LIGER, scVI
  new.reduction = "integrated.harmony",
  ndims        = 30
)

## Score the resulting embedding
obj <- IntegrationScore(
  obj,
  reduction = "integrated.harmony",
  ref.batch = obj@misc$reference,   # or any batch label
  batch     = "batch",
  K         = 10,
  n.cores   = 10
)

obj@misc$integration_effects[["integrated.harmony"]]
#>   BatchAlignment  CellIdentity      (per embedding dimension)
```

The two scores:

- **BatchAlignment** — higher = batches are better mixed (query batches match the reference's
  mixture structure).
- **CellIdentity** — higher = more genuine biological (multi-modal) structure is retained.

Their sum, **IntegrationScore = BatchAlignment + CellIdentity**, is the quantity the search
maximizes.

### Marginal vs. joint scoring

| Function                  | Model                                  | When to use                                  |
| ------------------------- | -------------------------------------- | -------------------------------------------- |
| `IntegrationScore`        | 1-D Gaussian mixture per dimension, aggregated by median | Default; robust and fast.   |
| `IntegrationScore.Joint`  | Single multivariate Gaussian mixture over all dimensions | Captures cross-dimension structure; set higher `K` (default 20). |

---

## Step 4 — Parameter search

`IntegrateRigor.ParameterS` automates the run-and-score loop across a grid of hyperparameters and
selects the combination with the highest IntegrationScore.

```r
## Build the grid — argument names must match the backend's parameter names.
param.df <- make.parameter.df(
  theta      = c(1, 2, 4),
  lambda     = c(1, 5)
)

obj <- IntegrateRigor.ParameterS(
  obj,
  parameter.df           = param.df,
  method                 = Harmony,
  batch                  = "batch",
  ndims                  = 30,        # embedding dimensions to compute
  ndims.score            = 30,        # dimensions to score over
  ref.batch              = NULL,       # auto-select if NULL
  K                      = 10,
  n.cores                = 5,
  use.batch.stable.genes = TRUE,       # adds "bsg" to the output reduction name
  subsample              = NULL,        # e.g. 0.2 to search on a subsample, then
  run.optimal            = TRUE         #   refit the winner on the full data
)
```

Results are stored on the object:

```r
obj@misc$parameter.search[["integrated.bsg.optimal.harmony"]]
#> Parameter | Reduction | BatchAlignment | CellIdentity | IntegrationScore

obj@misc$optimal.reduc        # which parameter set won
Embeddings(obj, "integrated.bsg.optimal.harmony")   # the winning embedding, ready for UMAP/clustering
```

For the joint multivariate score, use the drop-in `IntegrateRigor.ParameterS.Joint` (defaults to
`K = 20`, `n.cores = 1`).

### How `subsample` works

When `subsample` is set, the search runs on a random cell subset for speed. If `run.optimal = TRUE`,
the best parameter set is then re-fit on the **full** object so your final embedding uses all cells.