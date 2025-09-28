# _publish_modeling.R

library(bookdown)
library(gert)
library(glue)
library(fs)



repo_path <- "C:/Users/gajs/OneDrive - UW/UW/Courses/497 - Simulation/Textbook_Project/bookdown_497_simulation"
remote_expected <- "https://github.com/gajsivandran/Environmental_Modeling.git"
pages_url <- "https://gajsivandran.github.io/Environmental_Modeling/"

# 1) Sanity checks -----------------------------------------------------------
stopifnot(dir_exists(repo_path))
repo <- repo_path

# Ensure we’re in a git repo
git_info <- git_info(repo = repo)

# Ensure correct remote URL
remotes <- git_remote_list(repo = repo)
if (!any(remotes$name == "origin")) {
  git_remote_add(name = "origin", url = remote_expected, repo = repo)
} else {
  cur_url <- remotes$url[remotes$name == "origin"][1]
  if (!identical(cur_url, remote_expected)) {
    git_remote_set_url(remote = "origin", url = remote_expected, repo = repo)
  }
}

# Ensure we’re on 'main' (create/switch if needed)
branch <- tryCatch(git_branch(repo = repo)$name, error = function(e) NA_character_)
if (is.na(branch) || !nzchar(branch)) branch <- "main"
if (branch != "main") {
  # create local main if missing, then checkout
  branches <- git_branch_list(repo = repo)$name
  if (!"main" %in% branches) git_branch_create("main", checkout = FALSE, repo = repo)
  git_branch_checkout("main", repo = repo)
}

# 2) Build the book to docs/ ------------------------------------------------
old_wd <- getwd()
setwd(repo)
on.exit(setwd(old_wd), add = TRUE)

# Make sure docs exists and is not ignored
dir_create("docs")

# Optional: guard against docs being ignored
ignored <- tryCatch(git_check_ignore(paths = "docs/index.html", repo = repo), error = function(e) character())
if (length(ignored) > 0) {
  message("⚠️  It looks like docs/ is ignored by .gitignore. Remove that rule and re-run.")
  quit(status = 1)
}

# Build
render_book("index.Rmd", output_format = "bookdown::gitbook", output_dir = "docs")

# 3) Stage, commit, push -----------------------------------------------------
# Stage adds, modifications, and deletions inside docs/
git_add("docs", repo = repo)
# Handle deletions that git_add won't catch directly
# (gert stages deletions when you add the directory, but we can be explicit)
status <- git_status(repo = repo)

if (nrow(status) == 0) {
  message("ℹ️  No changes to publish. The site is already up to date.")
} else {
  n_add <- sum(status$status %in% c("new","modified","typechange"))
  n_del <- sum(status$status %in% c("deleted","renamed"))
  msg <- glue("Publish site: {n_add} changes, {n_del} deletions in docs/")
  git_commit(message = msg, repo = repo)
  # Pull (rebase) just in case, then push
  try(git_pull(repo = repo, fast_forward = TRUE), silent = TRUE)
  git_push(repo = repo)
  message("✅ Pushed changes.")
}

# 4) Final message -----------------------------------------------------------
cat(glue("\n✅ Published Modeling book to {pages_url} (branch main /docs)\n"))
cat("If the page looks stale, clear your browser cache or open a private window.\n")
