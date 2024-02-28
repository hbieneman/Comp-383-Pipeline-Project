import os
first_name = "Hailey"
last_name = "Bieneman"


output_directory = f"PipelineProject_{first_name}_{last_name}" # piece together name
os.system(f"mkdir {output_directory}") # making directory for output
os.chdir(output_directory)
