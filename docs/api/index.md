# API Reference

DeePhys is organised around a small set of value classes and handle classes that form a clear data flow:

```
SpikeData  →  RecordingProcessor  →  FeatureStore
                                          ↓
                              Classifier / Regressor / DimReducer
                                          ↓
                              ClassificationResult / RegressionResult / DimReductionResult
```

## Core classes

| Class | Type | Role |
|---|---|---|
| [SpikeData](SpikeData.md) | value | Raw Kilosort output container |
| [UnitData](UnitData.md) | value | Per-unit waveform and spike times |
| [RecordingProcessor](RecordingProcessor.md) | handle | Lazy computation engine for one recording |
| [FeatureStore](FeatureStore.md) | handle | Table-centric multi-recording feature container |
| [MLPipeline](MLPipeline.md) | static | ML helpers: train, CV, label pooling |

## Analysis classes

| Class | Role |
|---|---|
| [Classifier](Classifier.md) | Random forest / SVM / CNB / KNN classification |
| [Regressor](Regressor.md) | Random forest / SVM / KNN regression |
| [DimReducer](DimReducer.md) | UMAP / PCA / t-SNE dimensionality reduction |
| `CellTypeClassifier` | Transductive E/I labelling pipeline |
| `EIAnalyzer` | E/I activity traces and burst statistics |
| [Experiment](Experiment.md) | Orchestrates multi-recording analyses |

## Result classes

| Class | Role |
|---|---|
| [ClassificationResult](ClassificationResult.md) | Accuracy, confusion matrix, feature importances |
| [RegressionResult](RegressionResult.md) | R², RMSE, feature importances |
| `DimReductionResult` | Embedding coordinates and metadata |
| `MLResult` | Base class for ClassificationResult and RegressionResult (not instantiated directly) |

## Utility classes

| Class | Role |
|---|---|
| `FeatureAssembly` | Assembles feature vectors from raw unit/network data |
| [FeatureCatalog](FeatureCatalog.md) | Registry of named feature sets (`"full"`, `"core"`, custom) |
| [NormalizationPipeline](NormalizationPipeline.md) | Multi-step normalization (z-score / ComBat / clip / scale) with fit/transform pattern |
| `MLUtils` | Feature preprocessing helpers (NaN imputation, normalization, column pruning) |
| `StatisticalTest` | Group comparisons, ANOVA, paired tests, mixed-effects models, FDR correction |
| [GroupedCV](GroupedCV.md) | Hierarchy-aware grouped k-fold cross-validation |
| `ParentSpikeLoader` | Session-cached spike loading from parent (concatenated) recordings |
| `CacheManager` | On-disk computation caching |
| `AnalysisLog` | Structured logging for long batch runs |
| `RecordingDatabase` | SQLite-backed recording registry |
| `TrainedModel` | Serializable trained model wrapper |
