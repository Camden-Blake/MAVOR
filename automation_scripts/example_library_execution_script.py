from pathlib import Path
import os

mavor_exe = Path("..","build","mavor")
endf_loc = Path("..","build","test_files","ENDF8.0")
working_dir = Path("working_dir")
output_loc = Path("otf_files")

cmd = f"python process_endf_library.py {str(mavor_exe)} {str(endf_loc)} -d {str(working_dir)} -o {str(output_loc)}"

os.system(cmd)