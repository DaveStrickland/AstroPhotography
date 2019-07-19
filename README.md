# AstroPhotography

Python workflow for processing digital camera astrophotography


## Minimum Requirements

- Python 3.5+


## Optional Requirements

.. _pytest: http://pytest.org
.. _Sphinx: http://sphinx-doc.org

- `pytest`_ (for running the test suite)
- `Sphinx`_ (for generating documentation)


## Basic Setup

Install for the current user:

```bash
$ python -m pip install . --user
```

Run the application:

```bash
    $ python -m AstroPhotography --help
```

Run the test suite:

```bash
    $ pytest test/
```

Build documentation:

```bash
    $ sphinx-build -b html doc doc/_build/html
```
