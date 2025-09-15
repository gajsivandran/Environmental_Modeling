repo <- gert::git_find(".")
setwd(repo)

# Build
bookdown::clean_book(TRUE)
bookdown::render_book("index.Rmd","bookdown::gitbook")
if (!file.exists("docs/.nojekyll")) file.create("docs/.nojekyll")

# Stage + commit (ignore "nothing to commit")
gert::git_add(repo = repo, files = ".")
invisible(try(gert::git_commit(repo = repo,
                               message = sprintf("Publish: rebuild @ %s", Sys.time())), silent = TRUE))

# Detect current branch
br  <- gert::git_branch_list(repo = repo)
cur <- br$name[br$head]

# Push; if remote is ahead, fetch+merge then push again
tryCatch(
  gert::git_push(repo = repo, set_upstream = TRUE),
  error = function(e){
    if (grepl("contains commits that are not present locally", conditionMessage(e))) {
      message("Remote ahead; fetching and merging origin/", cur, " ...")
      gert::git_fetch(repo = repo, remote = "origin")
      gert::git_merge(repo = repo, commit = paste0("origin/", cur))
      gert::git_push(repo = repo, set_upstream = TRUE)
    } else stop(e)
  }
)
