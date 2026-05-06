# DimReducer

Static dimensionality reduction methods operating on plain feature tables. **No instantiation needed** — all methods are static.

## Example

```matlab
% UMAP on unit-level features
[X, unit_ids] = fs.unitMatrix('all');

opts.Method            = 'UMAP';
opts.NNeighbors        = 'auto';   % resolved to round(sqrt(N)), clamped to [5, 200]
opts.NormalizationGroups = fs.UnitTable.PlatingDate;
result = DimReducer.reduce(X, opts);

% Access embedding
coords = result.Embedding;   % (N x 2) 2D coordinates
scatter(coords(:,1), coords(:,2));

% PCA
result_pca = DimReducer.reduce(X, struct('Method', 'PCA', 'NDims', 10));

% t-SNE
result_tsne = DimReducer.reduce(X, struct('Method', 'tSNE', 'Perplexity', 30));
```

## Static methods

### `result = DimReducer.reduce(X, opts)`

Dimensionality reduction on a feature table.

| Argument | Type | Description |
|---|---|---|
| `X` | table | (N × F) feature table from `fs.unitMatrix()` etc. |
| `opts` | struct | Options (all fields optional; see below) |

Returns a `DimReductionResult` object.

## Options (`opts` struct fields)

All fields are optional.

| Field | Default | Description |
|---|---|---|
| `Method` | `"UMAP"` | `"UMAP"`, `"PCA"`, or `"tSNE"` |
| `NDims` | `2` | Number of output dimensions |
| `NNeighbors` | `"auto"` | UMAP k-NN neighbors; `"auto"` = `round(sqrt(N))`, clamped to [5, 200] |
| `MinDist` | `0.1` | UMAP minimum distance in embedding |
| `Spread` | `1.0` | UMAP effective scale of embedded points |
| `Perplexity` | `"auto"` | t-SNE perplexity; `"auto"` = `max(5, min(50, round(N/3)))` |
| `NormalizationGroups` | `[]` | (N × 1) group IDs for per-group z-score before reduction |
| `ColorFilePath` | `""` | Path to UMAP `colorsByName.properties` file (optional) |
| `Seed` | `[]` | RNG seed for reproducibility |

Auto parameters (`NNeighbors = "auto"`, `Perplexity = "auto"`) are resolved at runtime after data size N is known.

## `DimReductionResult` properties

| Property | Description |
|---|---|
| `Embedding` | (N × NDims) embedding coordinates |
| `Graph` | UMAP k-NN graph (sparse matrix); `[]` for PCA/tSNE |
| `Method` | String: `"UMAP"`, `"PCA"`, or `"tSNE"` |
| `Parameters` | Struct of all resolved options used |
| `Features` | Feature names used (string array) |
| `ElapsedTime` | Wall-clock time (seconds) |

## Notes

- UMAP requires the bundled `run_umap` function in `Toolboxes/umap/`.
- PCA and t-SNE use MATLAB built-ins (`pca`, `tsne`) from the Statistics and Machine Learning Toolbox.
- NaN values are imputed with the column mean before reduction (no train/test split — full dataset normalization).
- For joint embedding of DeePhys + external data (transfer learning), see [Tutorial 6](../tutorials/06_transfer_learning.md).
