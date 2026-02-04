


for f in *.pdf; do
  convert -density 300 "$f" -quality 90 "${f%.pdf}.png"
done



