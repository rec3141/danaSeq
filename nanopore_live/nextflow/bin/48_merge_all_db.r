library(DBI)
library(duckdb)
library(fs)
library(purrr)
library(stringr)
library(glue)

# where the hierarchy lives
root_dir <- "/data/project_QEI2025/out_dana/"

# --- discover all project / flowcell / barcode DBs ----------
dirs <- dir_ls(root_dir, recurse = F, type = "directory")
paths <- unname(unlist(sapply(dirs, function(x) dir_ls(x, recurse = 5, glob = "*.duckdb"))))

# parse components out of the path
meta <- tibble::tibble(
  path     = paths,
  project  = path_file(path_dir(path_dir(path_dir(path_dir(path))))),  # 3 levels up
  flowcell = path_file(path_dir(path_dir(path))),
  barcode  = path_file(path_dir(path))
)

# --- create / overwrite the merged DB ------------------------
dest <- "edna.duckdb"
if (file_exists(dest)) file_delete(dest)

con_dest <- dbConnect(duckdb(), dest, read_only = FALSE)

# --- merge all source DBs ------------------------------------
walk2(meta$path, seq_along(meta$path), function(srcdb, idx) {
  message("Merging ", srcdb)
  con_src <- dbConnect(duckdb(), srcdb, read_only = TRUE)
  
  schema_name <- glue("db_{idx}")  # unique schema name per attachment
  dbExecute(con_dest, glue("ATTACH '{srcdb}' AS {schema_name}"))
  
  for (tbl in dbListTables(con_src)) {
    message("  Processing table: ", tbl)
    
    # Format metadata tags for each row
    project  <- meta$project[idx]
    flowcell <- meta$flowcell[idx]
    barcode  <- meta$barcode[idx]
    
    if (idx == 1) {
      dbExecute(con_dest, glue(
        "CREATE TABLE {tbl} AS
       SELECT *, '{project}' AS project,
                 '{flowcell}' AS flowcell_src,
                 '{barcode}' AS barcode_src
       FROM {schema_name}.{tbl}"
      ))
    } else {
      dbExecute(con_dest, glue(
        "INSERT INTO {tbl}
       SELECT *, '{project}' AS project,
                 '{flowcell}' AS flowcell_src,
                 '{barcode}' AS barcode_src
       FROM {schema_name}.{tbl}"
      ))
    }
  }
  
  # Detach when done
  dbExecute(con_dest, glue("DETACH {schema_name}"))
})

dbDisconnect(con_dest, shutdown = TRUE)
message("Merged all data into: ", dest)
