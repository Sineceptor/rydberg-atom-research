# Manuscript draft

- **`DRAFT.md`**: the master draft. It opens with notes for the authors (delete before sharing), then the manuscript. Hidden `<!-- source: ... -->` comments give the data file behind every number.
- **`DRAFT.docx`**: Word export of `DRAFT.md` with native equations and embedded figures, for import into Google Docs. The source comments are dropped in this version. Rebuild it from the repository root:

  ```sh
  pandoc paper/DRAFT.md -f markdown+tex_math_dollars+raw_html -t docx --resource-path=. -o paper/DRAFT.docx
  ```

- **`figures/`**: publication figures regenerated from saved data by `make_figures.py`; the manifest `figures/FIGURES.md` lists sources and parameters.
- **`sections/`**: the per-section drafts from which `DRAFT.md` was assembled on 27 September. They are superseded by `DRAFT.md`, and some of their status lines are out of date. Edit `DRAFT.md`, not these.
