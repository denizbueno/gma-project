add_object_to_rda <- function(obj, rda_file, overwrite = FALSE) {
  if (!file.exists(rda_file)) {
    stop("The specified .RData file does not exist.")
  }

  # Create environments
  old_e <- new.env()
  new_e <- new.env()

  # Load existing objects into old_e
  load(file = rda_file, envir = old_e)

  # Get the name of the new object
  name_obj <- deparse(substitute(obj))

  # Assign the new object to new_e
  new_e[[name_obj]] <- obj

  if (overwrite) {
    # Overwrite existing objects with the same name
    invisible(sapply(ls(new_e), function(x) {
      assign(x, get(x, envir = new_e), envir = old_e)
    }))
    save(list = ls(old_e), file = rda_file, envir = old_e)
  } else {
    # Add new objects without overwriting
    invisible(sapply(ls(old_e), function(x) {
      assign(x, get(x, envir = old_e), envir = new_e)
    }))
    save(list = ls(new_e), file = rda_file, envir = new_e)
  }
}
