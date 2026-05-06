# FeatureCatalog

Central registry of named feature sets. **No instantiation needed** — all methods are static. Used by `FeatureStore` to filter feature columns at retrieval time.

## Built-in sets

| Set | Description |
|---|---|
| `"full"` | All features in each group (default, backwards-compatible) |
| `"core"` | Curated non-redundant subset: no algebraically-determined duplicates, interpretable neuroscientific meaning, low sensitivity to edge cases |

## Example

```matlab
% List available sets
sets = FeatureCatalog.availableSets();   % → ["core", "full"]

% Get feature names for a specific group and set
names = FeatureCatalog.features("ActivityFeatures", "core");

% Get all core features across all groups
all_core = FeatureCatalog.features("all", "core");

% Check if a feature name is known
tf = FeatureCatalog.isKnownFeature("FiringRate");   % → true

% Register a custom feature set for this session
FeatureCatalog.registerSet("pub", ["FiringRate", "T2Pdelay", "HalfWidth", "Synchrony"]);

% Use it in FeatureStore
[X, ids] = fs.unitMatrix('all', [], FeatureSet="pub");

% Remove when done
FeatureCatalog.removeSet("pub");
```

## Static methods

### `names = FeatureCatalog.features(group, set)`

Return feature names for a (group, set) combination.

| Argument | Description |
|---|---|
| `group` | Feature group name (e.g. `"ActivityFeatures"`, `"WaveformFeatures"`) or `"all"` for the union |
| `set` | `"full"` (default), `"core"`, or a custom registered set name |

For `set = "full"` and `group = "all"`, returns `string.empty` as a sentinel meaning "no column filter."

### `sets = FeatureCatalog.availableSets()`

Return all available set names as a string array, including custom registered sets.

### `tf = FeatureCatalog.isKnownFeature(name)`

Return `true` if `name` matches a known feature across any group and the `"full"` set.

### `FeatureCatalog.registerSet(name, feature_names)`

Register a custom named set for the current MATLAB session.

| Argument | Description |
|---|---|
| `name` | Valid MATLAB identifier string (cannot be `"core"` or `"full"`) |
| `feature_names` | String array of feature names to include |

Custom sets are cleared when the class is cleared (`clear FeatureCatalog`).

### `FeatureCatalog.removeSet(name)`

Remove a previously registered custom set. Built-in sets (`"core"`, `"full"`) cannot be removed.

## Notes

- Feature group names correspond to the column prefixes in `FeatureStore.UnitTable` (e.g. `"ActivityFeatures"`, `"WaveformFeatures"`, `"RegularityFeatures"`, `"ACG"`, `"FullACG"`, `"ReferenceWaveform"`, `"Catch22"`, `"NetworkFeatures"`).
- The `FeatureSet` option in `FeatureStore.unitMatrix()` and `Experiment` opts delegates to `FeatureCatalog.features()` internally.
- Using `"core"` is recommended for publication-quality results: it removes redundant and unstable features.
