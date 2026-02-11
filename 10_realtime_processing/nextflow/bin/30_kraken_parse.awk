# kraken_parse.awk
BEGIN {
  FS = "\t"
}
{
  id = $2
  if ($1 == "C") {
    tax = $3
  gsub(/ /, "_", tax)
  print id "\t" tax
  } else {
    tax = "Unknown"
  }
}
