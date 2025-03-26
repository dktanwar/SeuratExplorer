# Before build R package for CRAN

## 1. DESCRIPTION file

Remove line bellow:

```
Remotes:
    bioc::release/ComplexHeatmap,
    github::immunogenomics/presto,
    bioc::release/MAST,
    bioc::release/limma,
    bioc::release/DESeq2
# may be also need to remove Suggests part for non-CRAN packages.
```

Hopefully add doi info to DESCRIPTION.

## 2. Remove files

Remove files bellow:

```
# 1. this file: before-submit-to-CRAN.md
# 2. CITATION.cff
```

## 3. Build the package

Build - More- Build Source Package.



Attention： update the CRAN package each two month.