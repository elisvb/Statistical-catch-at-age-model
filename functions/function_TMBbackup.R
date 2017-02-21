##### TMBbackup-----
# read a cpp file, without running it (for the backup function)
read_cpp_file <- function(path) {
  paste0(paste0(readLines(path), collapse = "\n"), "\n")
}

# make a backup of a cpp file (for when model is actually working)
TMBbackup = function(path,filename){
  input=read_cpp_file(path)
  write(input, file = filename)
}