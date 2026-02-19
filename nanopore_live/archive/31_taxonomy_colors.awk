# taxonomy_colors.awk

function ord(c,    chars) {
  chars = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
  return index(chars, c)
}

function hash_color(tax,  i, ch, h) {
  h = 0
  for (i = 1; i <= length(tax); i++) {
    ch = substr(tax, i, 1)
    h += ord(ch) * i
  }
  return sprintf("#%06x", h % 16777215)
}

BEGIN {
  OFS = ","
  print "Name", "Taxonomy", "Color"
}

{
  name = $1
  tax = $2
  if (!(tax in color)) {
    color[tax] = hash_color(tax)
  }
  print name, tax, color[tax]
}
