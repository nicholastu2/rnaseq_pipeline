help(
[[
This module loads the Brent lab rnaseq_pipeline.
]])

local pkgName = myModuleFullName()
local fullVersion = myModuleVersion ()
local base = "/opt/apps/labs/mblab/software/rnaseq_pipeline/1.0"

-- display this info with module whatis
whatis("Name: " ..pkgName)
whatis("Version: " ..fullVersion)
whatis("Description: The Brent lab RNA-Seq pipeline")
whatis("Location: " ..base)
whatis("URL: https://gitlab.com/brentlab/rnaseq_pipe")

-- load main tools (see system requirements on the gitlab wiki)
depends_on("r/3.5.3-python-3.6.5")
depends_on("igv/2.4.7")
depends_on("py-pandas/0.24.1-python-3.6.5")
depends_on(" py-numpy/1.16.2-python-3.6.5")
depends_on("py-pyyaml/3.13-python-3.6.5")
depends_on("py-openpyxl/2.4.5-python-3.6.5")

-- message to print on module load
if(mode() == 'load') then
   LmodMessage("\nBrent lab rnaseq_pipeline is loaded!\n\tPlease see https://gitlab.com/brentlab/rnaseq_pipe/-/wikis/home for usage instructions.\n")
end

-- message to print on module unload
if (mode() == "unload") then
   LmodMessage("\nmodule rnaseq_pipeline and all dependencies automatically loaded by rnaseq_pipeline are unloaded.\n")
end

-- set environmental variables for R and python dependencies
setenv("R_LIBS_USER", "/opt/apps/labs/mblab/software/rnaseq_pipeline/v1.0_r_dependencies")
setenv("R_INCLUDE_DIR", "/opt/apps/labs/mblab/software/rnaseq_pipeline/v1.0_r_dependencies")
--setenv("PYTHONPATH", "/opt/apps/labs/mblab/software/rnaseq_pipeline/v1.0_python_dependencies")

-- prepend the following to $PATH
prepend_path("PATH", base)
prepend_path("PATH", pathJoin(base, "tools"))








