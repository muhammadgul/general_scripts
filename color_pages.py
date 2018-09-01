import numpy as np

# gs -o - -sDEVICE=inkcov /path/to/your.pdf &> tmp.txt

pages = ""
counter = 0
colored = 0
with open("tmp.txt",'r') as f:
        for line in f:
                line = line.strip()
                
                if ("CMYK OK" in line):
                        counter += 1
                        if not "0.00000  0.00000  0.00000  0." in line:
                                pages += "," + str(counter)
                                colored+=1
print(colored)
print(pages)
