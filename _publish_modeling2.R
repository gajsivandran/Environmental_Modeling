# _publish_calculus.R
# Bookdown -> GitHub Pages (docs/) with auto-fixes for the Calculus book

EXPECTED_REMOTE <- "https://github.com/gajsivandran/Environmental_Modeling.git"
PAGES_URL       <- "https://gajsivandran.github.io/Environmental_Modeling/"
BRANCH_OVERRIDE <- Sys.getenv("BRANCH_OVERRIDE", unset = NA_character_)  # e.g., "main2"

suppressPackageStartupMessages({
  if (!requireNamespace("gert", quietly = TRUE)) install.packages("gert")
  if (!requireNamespace("bookdown", quietly = TRUE)) install.packages("bookdown")
  if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml")
})

# ---------- helpers ----------
nuke_knitr_caches <- function(root = ".") {
  # remove common knitr cache dirs (e.g., bookdown-demo_cache, cache, *_cache)
  pats <- c("*_cache", "cache")
  killed <- character(0)
  for (p in pats) {
    hits <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    hits <- hits[basename(hits) == p | grepl("_cache$", basename(hits))]
    for (h in unique(hits)) {
      if (dir.exists(h)) {
        unlink(h, recursive = TRUE, force = TRUE)
        killed <- c(killed, h)
      }
    }
  }
  if (length(killed)) message("Removed knitr caches: ", paste(killed, collapse = ", "))
  invisible(killed)
}


ensure_bookdown_outputs_to_docs <- function() {
  yml_path <- "_bookdown.yml"
  if (!file.exists(yml_path)) {
    yaml::write_yaml(list(output_dir = "docs"), yml_path)
    message("Created _bookdown.yml with output_dir: 'docs'")
    return(invisible(TRUE))
  }
  y <- tryCatch(yaml::read_yaml(yml_path), error = function(e) NULL)
  if (is.null(y)) stop("Could not read _bookdown.yml. Fix its syntax and re-run.")
  if (is.null(y$output_dir) || !identical(y$output_dir, "docs")) {
    y$output_dir <- "docs"
    yaml::write_yaml(y, yml_path)
    message("Set output_dir: 'docs' in _bookdown.yml")
  }
  invisible(TRUE)
}

ensure_nojekyll <- function() {
  if (!dir.exists("docs")) dir.create("docs", recursive = TRUE, showWarnings = FALSE)
  if (!file.exists("docs/.nojekyll")) {
    file.create("docs/.nojekyll")
    message("Created docs/.nojekyll")
  }
}

ensure_gitignore_allows_docs <- function() {
  gi <- ".gitignore"
  if (!file.exists(gi)) return(invisible(TRUE))
  lines <- readLines(gi, warn = FALSE)
  # remove any line that is exactly 'docs' or '/docs'
  bad <- grepl("^\\s*(/?)docs/?\\s*$", lines)
  if (any(bad)) {
    lines <- lines[!bad]
    writeLines(lines, gi)
    message("Removed 'docs' from .gitignore")
  }
  invisible(TRUE)
}

detect_branch <- function(repo) {
  # 1) explicit override wins
  if (!is.na(BRANCH_OVERRIDE) && nzchar(BRANCH_OVERRIDE)) return(BRANCH_OVERRIDE)
  # 2) try gert::git_info()
  cur <- tryCatch(gert::git_info(repo)$branch, error = function(e) NA_character_)
  # Guard against NULL/length 0/NA/""
  if (is.null(cur) || length(cur) != 1L || is.na(cur) || identical(cur, "")) {
    return("main")
  }
  cur
}

push_with_merge_if_needed <- function(repo, branch) {
  push_once <- function() {
    gert::git_push(
      repo = repo, remote = "origin",
      refspec = paste0("HEAD:refs/heads/", branch),
      set_upstream = TRUE
    )
  }
  tryCatch(
    { push_once(); TRUE },
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("not present locally|non-fast-forward|contains commits", msg)) {
        message("Remote ahead; fetching and merging origin/", branch, " ...")
        gert::git_fetch(repo = repo, remote = "origin")
        gert::git_merge(repo = repo, commit = paste0("origin/", branch))
        st <- gert::git_status(repo = repo)
        if (any(st$status == "conflicted")) {
          stop("Merge conflicts detected. Resolve them, commit, then re-run.")
        }
        push_once(); TRUE
      } else stop(e)
    }
  )
}

# ---------- main ----------
repo <- tryCatch(gert::git_find("."), error = function(e)
  stop("Not inside a git repo. Open the CALCULUS book project and run again.")
)
setwd(repo)
if (!file.exists("index.Rmd")) stop("index.Rmd not found in: ", repo)

# Ensure origin -> Calculus repo (auto-fix)
rem <- gert::git_remote_list(repo = repo)
if (nrow(rem) == 0) {
  gert::git_remote_add(remote = "origin", url = EXPECTED_REMOTE, repo = repo)
  message("Added origin: ", EXPECTED_REMOTE)
} else {
  origin_url <- rem$url[rem$name == "origin"][1]
  if (!identical(origin_url, EXPECTED_REMOTE)) {
    gert::git_remote_set_url(remote = "origin", url = EXPECTED_REMOTE, repo = repo)
    message("Updated origin to: ", EXPECTED_REMOTE)
  }
}

# Output dir/docs housekeeping
ensure_bookdown_outputs_to_docs()
ensure_nojekyll()
ensure_gitignore_allows_docs()

# Branch
branch <- detect_branch(repo)
message("Repo   : ", repo)
message("Origin : ", EXPECTED_REMOTE)
message("Branch : ", branch)

# Build
bookdown::clean_book(TRUE)
nuke_knitr_caches(".")   # <-- add this line
bookdown::render_book("index.Rmd", "bookdown::gitbook")
ensure_nojekyll()
if (!file.exists("docs/index.html")) stop("Build failed: docs/index.html not found")

# Commit + push (auto-handle remote-ahead)
gert::git_add(repo = repo, files = ".")
invisible(try(
  gert::git_commit(
    repo = repo,
    message = sprintf("Publish (Calculus): rebuild @ %s",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  silent = TRUE
))
push_with_merge_if_needed(repo, branch)

message("✅ Published Calculus book to ", PAGES_URL, " (branch ", branch, ")")
message("If the page looks stale, hard refresh the browser.")
