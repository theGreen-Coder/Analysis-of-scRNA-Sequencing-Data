import os
import sys

sysArgs = sys.argv
sysArgs.remove("run.py")
args = " ".join(sysArgs)

os.system(("python3 main.py "+args))
os.system("python3 trajectory.py "+args)