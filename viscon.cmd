@echo off

set CURRDIR=%~dp0
set PYTHONDIR=%CURRDIR%\WinPython
set PYTHONBIN=%PYTHONDIR%\python-3.6.1.amd64
set GMSHBIN=%PYTHONDIR%\tools

set PYTHONPATH=%PYTHONBIN%
set PATH=%PYTHONBIN%;%GMSHBIN%;%PATH%

python viscon.py