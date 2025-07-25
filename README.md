# GeoFetcheroo

GeoFetcheroo is a playful yet powerful command‑line utility that bridges GEO and SRA to give you all the sequencing run details you need.

## Additional Details

- **One‑step workflow**: From a GEO Series (`GSE`) to SRR runs and human‑readable sample names in a single script.  
- **No manual scraping**: Automates FTP, HTML or CSV parsing by using R’s **GEOquery** and NCBI’s **EDirect** tools.  
- **Tab‑delimited output**: Easy to import into spreadsheets, R/Python data frames or downstream pipelines.

## Dependencies

- **R** (≥ 4.0)  
- **GEOquery** R package  
  ```r
  install.packages("GEOquery")
```
