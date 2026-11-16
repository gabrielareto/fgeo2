# Example datasets

## Barro Colorado Island environmental extraction

The environmental extraction tutorial uses the BCI 50-ha forest dynamics plot,
with local coordinates spanning 1,000 m west to east and 500 m south to north.

| File | Contents |
| --- | --- |
| `census1.txt` | Existing CTFS tab-separated stem records for census 1, nominally 1982. There are 263,934 rows, with measurement dates from March 1981 to July 1983. |
| `census2.txt` | Existing CTFS tab-separated stem records for census 2 in 1985. There are 299,447 rows. |
| `bci-elevation.csv` | Complete historical BCI elevation grid with 20,301 rows and columns `x`, `y`, `z`. Coordinates and elevation are in metres. |
| `bci-soils-wolf-2015.csv` | Complete Wolf et al. soil dataset with 625 sampling locations and all 30 original columns. |

### Tree censuses

The two census files are the historical CTFS text exports already present in
this repository. Their original release version is undocumented here. The
curated BCI archive describes the study, census history, coordinate conventions
and data dictionaries. Credit the BCI investigators Richard Condit, Rolando
Pérez, Salomón Aguilar, Suzanne Lao, Robin Foster and Stephen P. Hubbell.
[BCI data archive, Condit et al. 2019](https://doi.org/10.15146/5xcp-0d46).

In these exports, `PX` and `PY` are plot coordinates in metres, `Mnemonic` is
the species code, and `StemID` identifies a stem within a census. `TreeID` can
repeat because a tree can have several stems. Preserve `Tag`, `StemTag` and
`Quadrat` as character columns when reading the files, so leading zeros survive.
The pair `{Census, StemID}` uniquely identifies each supplied census record.

### Elevation

`bci-elevation.csv` was converted from `BCI elev.TXT` in Gabriel Arellano's
existing BCI analysis files, with a source-file modification date of 19 April
2017. The conversion changed the delimiter to comma and renamed `elev` to `z`.
All numerical values and rows were retained. The grid contains every combination
of `x = 0, 5, ..., 1000` and `y = 0, 5, ..., 500`.

The values were verified against the published `BCIelev.tsv` in the BCI archive
by Condit et al. Normalizing the historical table to tab-separated numeric text
reproduced the archive file's full SHA-256 checksum,
`0869ea08d2488b1e5c67d6e78201c6d0e938123be74ae190210d411aefd683b4`.
The archive defines `x` and `y` as distances from the west and south plot borders,
and elevation relative to sea level. Its methods refer to Harms et al. (2001),
*Habitat associations of trees and shrubs in a 50-ha neotropical forest plot*.
The 5 m spacing describes the supplied grid; its survey and preprocessing history
determine the underlying accuracy. The Dryad release is **CC0 1.0**.
[BCI data archive](https://doi.org/10.15146/5xcp-0d46) ·
[Published elevation file metadata](https://datadryad.org/api/v2/files/148947).

### Soil chemistry

Wolf, J. A., Hubbell, S. P., Fricker, G. A. and Turner, B. L. (2015),
*Geospatial observations on tropical forest surface soil chemistry*,
Ecology 96, 2313. [Paper DOI](https://doi.org/10.1890/15-0558.1).

This file is an unchanged copy of the deposited `soils.csv`, with a more
descriptive filename. It contains 2010 surface-soil measurements, local plot
coordinates `gx`, `gy`, UTM coordinates, nutrients, carbon, nitrogen and pH.
Keep the extraction methods in the variable names. For example, `M3Ca` is
Mehlich-3 calcium and `Ca` is BaCl2-extractable calcium. The header contains
a leading space in ` utm_n`; the tutorial removes header whitespace on import.
[Variable definitions and quality-control notes](https://esapubs.org/archive/ecol/E096/204/metadata.php).

The publisher's dataset is released under **CC0 1.0**. Its original metadata
also explicitly states that the data have no copyright restrictions. Retain
the source citation and identify any derived or interpolated products.
[Dataset and licence](https://doi.org/10.6084/m9.figshare.3562863) ·
[CC0 terms](https://creativecommons.org/publicdomain/zero/1.0/).

The unchanged file has MD5 `d4deb84fa3a2300097efc07a7a42700c`. Its original
carriage-return line endings are retained. [Original CSV download](https://ndownloader.figshare.com/files/5635584).

Dataset-specific terms apply independently of the licence for the R code.
