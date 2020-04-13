import sys
from pathlib import Path

path = Path('../pytrack/rjtd/track2019100900.txt')
f = path.open()
l = f.readline().split()
print(l)
