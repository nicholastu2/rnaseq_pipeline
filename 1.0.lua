help(
[[
This module loads the Brent lab rnaseq_pipeline.
]])

local pkgName = myModuleFullName()
local fullVersion = myModuleVersion ()
local base = "/opt/apps/labs/mblab/software/rnaseq_pipeline/1.0"
local tools = pathJoin(base, "tools")
local metadata_database = "/lts/mblab/Crypto/rnaseq_data/database-files" 

-- display this info with module whatis
whatis("Name: " ..pkgName)
whatis("Version: " ..fullVersion)
whatis("Description: The Brent lab RNA-Seq pipeline")
whatis("Location: " ..base)
whatis("URL: https://gitlab.com/brentlab/rnaseq_pipe")

-- load main tools (see system requirements on the gitlab wiki)
depends_on("miniconda")
depends_on("igv/2.4.7")
depends_on("samtools")

-- message to print on module load
if(mode() == 'load') then
   LmodMessage("\nBrent lab rnaseq_pipeline is loaded!\n\tPlease see https://github.com/BrentLab/rnaseq_pipeline/wiki for usage instructions.\n")
end

-- prepend the following to $PATH
prepend_path("PATH", base)
prepend_path("PATH", tools)

-- set environmental variables
setenv("METADATA", metadata_database)
setenv("PYTHONPATH", tools)


-- see test_conda_env in the rnaseq_pipeline code repo. This simply tests whether the directory $HOME/.conda/envs exists. If it does not, the rnaseq_pipeline_env directory is cloned by conda
-- into $HOME/.conda/envs.
-- see the miniconda module help for more details
if (not isDir("$HOME/.conda/envs/rnaseq_pipeline")) then
  execute{cmd="test_conda_env",modeA={"load"}}
end

-- load/unload virtual environment located in $HOME/.conda/envs/rnaseq_pipeline
execute{cmd="conda activate rnaseq_pipeline",modeA={"load"}}

-- message to print on module unload
if (mode() == "unload") then
   LmodMessage("\nmodule rnaseq_pipeline and all dependencies and environments loaded by rnaseq_pipeline are unloaded.\n")
end
