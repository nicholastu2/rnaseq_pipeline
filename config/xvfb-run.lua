help(
[[
This module loads xvfb-run. This module executes programs that require a monitor output on a headless system
]])

local pkgName = myModuleFullName()
local base = "/opt/apps/labs/mblab/software/xvfb"

-- display this info with module whatis
whatis("Name: " ..pkgName)
whatis("Description: this module was installed for to create igv snapshots on the cluster")
whatis("Location: " ..base)
whatis("URL: http://manpages.ubuntu.com/manpages/xenial/man1/xvfb-run.1.html")

-- message to print on module load
if(mode() == 'load') then
   LmodMessage("\nxvfb-run is loaded\n")
end

prepend_path("PATH", base)
