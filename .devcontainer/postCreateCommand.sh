# Activate the codespaces configuration for vscode
/usr/bin/ln -s /workspaces/GEOS/.vscode-codespaces /workspaces/GEOS/.vscode
# Activate the appropriate submodules
git submodule init
git submodule update
# Load the pretty printer for LvArray
echo "source /workspaces/GEOS/src/coreComponents/LvArray/scripts/gdb-printers-shallow.py" > ~/.gdbinit
