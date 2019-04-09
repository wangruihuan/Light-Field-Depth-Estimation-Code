@echo off
for %%i in (*.pfm) do ( python2 E:\DepthEstimation\tool\evaluation-toolkit-master\evaluation-toolkit-master\source\convert_pfm2png.py %%i parameters.cfg %%~ni.png)
cmd