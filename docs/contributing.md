# Contributing

Contributions are welcome — bug fixes, new features, additional tutorials, and documentation improvements.

## Reporting issues

Use the [GitHub issue tracker](https://github.com/hornauerp/DeePhys/issues). Please include:
- MATLAB version and OS
- Steps to reproduce
- Any error messages or stack traces

## Development setup

```bash
git clone https://github.com/hornauerp/DeePhys.git
cd DeePhys
```

Run `startup.m` in MATLAB from the repo root to add all paths.

## Running tests

```matlab
% From the repo root in MATLAB — standard MATLAB test runner:
cd Tests
runtests('.')

% Or use the project's own runner:
run_all_tests
```

Tests live in `Tests/` and follow the `test_*.m` naming convention.

## Pull requests

1. Fork the repository and create a feature branch
2. Make changes, add or update tests where applicable
3. Ensure all existing tests still pass
4. Open a pull request against `main` with a description of what changed and why

## Code style

- Class files go in `Classes/@ClassName/`
- Utility functions go in `Functions/`
- Each public method should have a one-line `%` help comment immediately after the function signature
- Value classes for data containers, handle classes for stateful computation engines

## Documentation

Documentation lives in `docs/` and is built with MkDocs + Material theme.

```bash
pip install mkdocs-material
mkdocs serve   # preview at http://127.0.0.1:8000
```

API pages are in `docs/api/`, one file per class. Follow the existing template (Example → Constructor → Properties → Methods).
