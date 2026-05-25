# Backend Note

The maintained storage backends are:

- `lance`
- `zarr`

The maintained corpus topologies are:

- `aggregate`
- `federated`

Recommended default:

- aggregate Lance

Supported secondary routes:

- federated Lance
- aggregate Zarr
- federated Zarr

Use aggregate Lance for the normal create, append, canonicalize, and load path.
Use Zarr when chunked array artifacts or node-local staging are operationally
useful. Use federated topology when whole-dataset isolation or recomposition is
more important than one shared aggregate matrix object.

Historical backend experiments are not part of the current mainline docs.
