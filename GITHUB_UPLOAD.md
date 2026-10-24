# Uploading this repository to GitHub

Target repository:

`https://github.com/ygwang2018/Sparsity-Aware-Estimation`

If the repository already exists, upload the **contents of this folder** to the repository root (not the enclosing folder itself). The root should contain `DESCRIPTION`, `NAMESPACE`, `README.md`, `R/`, `man/`, `inst/`, and `tests/`.

After upload, users can install with:

```r
remotes::install_github("ygwang2018/Sparsity-Aware-Estimation")
```

Recommended local validation before tagging a release:

```r
install.packages(c("devtools", "testthat"))
devtools::document()
devtools::test()
devtools::check()
```
