# AstroPhotography

Python workflow for processing digital camera astrophotography


## Minimum Requirements

See `requirements.txt` for full dependency list.
- Python 3.5+
- `numpy` and `matplotlib`
- `PyYAML`
- `rawpy` (python binding and interface to libraw)


## Optional Requirements

- `pytest` http://pytest.org (for running the test suite)
- `Sphinx` http://sphinx-doc.org (for generating documentation)


## Basic Setup

Install:

```bash
$ python3 -m pip install . --requirements=requirements.txt
```

Run the application:

```bash
    $ python3 -m AstroPhotography --help
    # OR
    $ python3 dksraw --help
```

Run the test suite:

```bash
    $ pytest test/
```

Build documentation:

```bash
    $ sphinx-build -b html doc doc/_build/html
```
