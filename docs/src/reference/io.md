# Data loading

Spectral reference tables bundled in `data/` are loaded through the
utilities below. The API is simple enough that downstream users can
drop their own provenance-commented CSVs into `data/` and pick them up
with the same loader.

```@docs
SpectralTable
linterp
data_path
load_spectral_csv
```
