# Sharing the ICI scorer with colleagues (for feedback)

There are two ways to let colleagues try the app. The **static bundle** needs
no R and no server account: it runs entirely in their browser.

## Option 1 — Static bundle (no install, recommended for feedback)

The app has been compiled with [shinylive](https://posit-dev.github.io/shinylive/)
into a self-contained static website (a folder of HTML/JS/WebAssembly). It runs
R *inside the browser* via webR, so colleagues need nothing but a modern browser
(Chrome, Edge, Firefox, Safari). First load fetches ~30-60 MB of R packages and
is slow (20-40 s); it is cached afterwards.

The built site is in **`ici_shinylive_site/`** (also zipped as
`ici_shinylive_site.zip`).

To give colleagues a link, host the folder on any static host:

- **Netlify Drop (no account needed):** go to <https://app.netlify.com/drop>
  and drag the `ici_shinylive_site` folder onto the page. You get a public URL
  in seconds; share it. (Free, temporary unless you claim it with an account.)
- **GitHub Pages:** copy the folder into a repo's `docs/` (or a `gh-pages`
  branch), enable Pages, and share the Pages URL.
- **Any web server / S3 bucket:** upload the folder; serve `index.html`.

To preview it locally before sharing:

```bash
cd ici_shinylive_site
python3 -m http.server 8000
# open http://localhost:8000
```

Opening `index.html` directly from disk (file://) will NOT work — it must be
served over http(s).

## Option 2 — shinyapps.io (a durable, hosted link)

If you want a permanent, robust link and are happy to make a free account:

1. Create an account at <https://www.shinyapps.io> and get your token
   (Account -> Tokens).
2. In R, from the project root:

   ```r
   install.packages("rsconnect")
   rsconnect::setAccountInfo(name="<you>", token="<token>", secret="<secret>")
   rsconnect::deployApp("app", appFiles = c(
     "app.R", "ici_locked_reference.R", "score_cohort.R",
     "references/ICI_v1.rds", "example_cohort.csv"))
   ```

   (Leave out `make_app_assets.R` - it needs Bioconductor packages only used to
   rebuild the reference, not to run the app.)

This runs the app on a real R server, so every feature works exactly as it does
locally.

## What colleagues can do

- **Patient report:** pick a patient (or upload a one-sample file), read the
  call, the 3-layer percentile gauge, and the cascade drill-down.
- **Cohort explorer:** load their own cohort (CSV / TSV / Excel, gene x sample)
  or the bundled example; see archetype proportions, layer distributions, the
  per-patient beeswarm, the table, and download a cohort PDF.
- Use **Load example cohort** to try it instantly with no file.

**Research use only.** The app is investigational and ships only public example
data (an IMvigor210 subset); it contains no patient data.
