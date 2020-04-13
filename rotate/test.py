import sys

s = sys.stdin.readline().strip("\n").split(" ")
print(s)
print(type(s))
p1 = s[0]
p2 = int(s[1])
p3 = float(s[2])
print(type(p1),type(p2),type(p3))
