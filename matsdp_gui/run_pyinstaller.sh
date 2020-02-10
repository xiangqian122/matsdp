#!/bin/bash
pyinstaller -F -w --distpath ../ -n matsdp_gui.exe --onefile --hidden-import="sklearn.neighbors.typedefs" --hidden-import="sklearn.utils._cython_blas"  --hidden-import="sklearn.neighbors.quad_tree" --hidden-import="sklearn.tree._utils" --hidden-import="PyQt5" gui_main.py
