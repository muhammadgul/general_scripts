```
# To find the color pages in a pdf file

gs -o - -sDEVICE=inkcov /path/to/your.pdf &> tmp.txt
python color_pages.py tmp.txt






```
