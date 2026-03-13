# Genetic Association Heatmap

An interactive R Shiny application for visualising the **beta values** of
genetic variant–study associations.  The heatmap is built with
[ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap/) and made
interactive with
[InteractiveComplexHeatmap](https://bioconductor.org/packages/InteractiveComplexHeatmap/).

Data can be sourced from two REST APIs:

| Source | URL |
|--------|-----|
| GWAS Catalog REST API v2 | <https://www.ebi.ac.uk/gwas/rest/api/v2> |
| HuGeAMP BioIndex REST API | <https://bioindex.hugeamp.org/api/bio> |

---

## Features

- **Variant search** — query either API by rsID to find all available studies
  and their beta values.
- **Bulk variant input** — paste a list of rsIDs (one per line) to add many
  variants in a single operation, with an option to limit results to studies
  already shown in the heatmap.
- **Flexible single-variant addition** — add a variant to every study currently
  shown in the heatmap, to every study returned by the search, or to
  individually selected studies.
- **Manual entry** — type a variant ID + study ID to add a single cell.
- **Auto-fill** — whenever a new variant or study is added, the app
  automatically fetches all missing beta values so the matrix is always
  complete (variants × studies).
- **Interactive heatmap** — click a cell to see variant, study, and beta
  details; brush to inspect a region.
- **PNG download** — export the current heatmap as a high-resolution PNG.
- **Data table** — inspect all raw beta values in a searchable table.

---

## Quick Start (RStudio)

### 1 — Prerequisites

Install R (≥ 4.1) and RStudio, then install the required packages once:

```r
# CRAN packages
install.packages(c("shiny", "httr", "jsonlite", "DT", "circlize"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ComplexHeatmap", "InteractiveComplexHeatmap"))
```

### 2 — Get the App

Clone the repository or download it as a ZIP:

```bash
git clone https://github.com/manning-lab/shiny_genetic_association_heatmap.git
```

### 3 — Run from RStudio

Open the project folder in RStudio, then in the R console run:

```r
shiny::runApp()
```

RStudio will detect `app.R` automatically and launch the app in your browser.

Alternatively, open `app.R` in the editor and click the **Run App** button in
the top-right corner of the editor pane.

### 4 — Run from any R session

```r
shiny::runApp("/path/to/shiny_genetic_association_heatmap")
```

---

## Walkthrough

1. **Select a Data Source** — choose *GWAS Catalog v2* or *HuGeAMP BioIndex*
   in the sidebar.

2. **Search for a single variant** — enter an rsID (e.g. `rs7903146`) and click
   **Search**.  The results table lists every study with a known association and
   its beta value.  Use one of the three **Add** buttons to add it to the
   heatmap.

3. **Add many variants at once (Bulk Variant Input)** — paste a list of rsIDs,
   one per line, into the *Bulk Variant Input* textarea and click
   **Add All Variants**.
   - Lines starting with `#` are treated as comments and ignored.
   - Duplicate IDs are automatically deduplicated.
   - Enable **Limit to studies already in heatmap** to fetch betas only for
     the studies currently represented, so you can expand the variant axis
     without adding new studies.

4. **Manual entry** — type both a Variant ID and a Study ID and click
   **Add to Heatmap** to add an individual cell.

5. **Auto-fill** — after any addition the app fetches missing beta values so
   the heatmap matrix is always rectangular (all variants × all studies).

6. **Explore the heatmap** — click a cell to read variant / study / beta
   details; click-drag to select a region and see summary statistics.

7. **Download** — click **Download PNG** to save the current heatmap.

---

## Project Structure

```
shiny_genetic_association_heatmap/
├── app.R                  # Main Shiny application (UI + server)
├── R/
│   ├── api_utils.R        # GWAS Catalog & HuGeAMP API helper functions
│   └── heatmap_utils.R    # ComplexHeatmap construction utilities
└── README.md
```

---

## API Notes

### GWAS Catalog REST API v2

Endpoint used: `GET /associations?filter=variantId:{rsId}&size=500`

The response is HAL-JSON (`_embedded.associations`).  Each association object
contains `beta`, `pValue`, `study.accessionId`, and `efoTraits`.

### HuGeAMP BioIndex REST API

Endpoint used: `GET /api/bio/query/associations?q={rsId}&limit=1000`

The response contains a `data` array.  Column names are detected dynamically
so the app handles minor API variations gracefully.

---

## License

MIT
