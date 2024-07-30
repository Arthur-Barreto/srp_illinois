from pathlib import Path
import subprocess
import numpy as np

# Define paths and parameters
MAIN_PATH = Path.cwd().parent
STAPL_PATH = MAIN_PATH.parent / "stapl" / "stapl-developer" / "examples" / "bls"
# /home/abarreto/arthur/stapl/stapl-developer/examples/bls
BUILD_FOLDER = MAIN_PATH / "build"
HOST_FILE_FOLDER = MAIN_PATH / "host_file"

print(STAPL_PATH)

execs = ["bls_omp", "bls_mpi", "bls_stapl"]
# execs = ["bls_mpi"]
num_nodes = [2**i for i in range(6)]
# n_points = 5000
n_points = np.arange(5000, 65000 + 1, 5000)
num_samples = 32

for exec in execs:

    EXECUTABLE = BUILD_FOLDER / exec

    if exec == "bls_stapl":
        EXECUTABLE = STAPL_PATH / exec

    for node in num_nodes:

        for n in n_points:

            # Check if executable exists
            if not EXECUTABLE.exists():
                print(f"Executable not found at {EXECUTABLE}")
            else:
                if exec == "bls_omp":
                    command = f"{EXECUTABLE} {n} {node} {num_samples}"
                elif exec == "bls_mpi":  # exec == "bls_mpi"
                    command = f"mpiexec -f {HOST_FILE_FOLDER} -n {node} {EXECUTABLE} {n} {num_samples}"

                else:
                    command = f"mpiexec -f {HOST_FILE_FOLDER} -n {node} {EXECUTABLE} {n} {num_samples}"

                print(f"Running command: {command}")

                try:
                    result = subprocess.run(
                        command,
                        shell=True,
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                    )
                    if exec == "bls_omp":
                        with open("log_bls_omp.txt", "a") as f:
                            f.write(f"{command} \n")
                            f.write(result.stdout.decode())
                            f.write(result.stderr.decode())
                            f.write("\n")
                    elif exec == "bls_mpi":
                        with open("log_bls_mpi.txt", "a") as f:
                            f.write(f"{command} \n")
                            f.write(result.stdout.decode())
                            f.write(result.stderr.decode())
                            f.write("\n")
                    else:
                        with open("log_bls_stpal.txt", "a") as f:
                            f.write(f"{command} \n")
                            f.write(result.stdout.decode())
                            f.write(result.stderr.decode())
                            f.write("\n")
                except subprocess.CalledProcessError as e:
                    print("Error running command:", e)
                    print("Command output:", e.stdout.decode())
                    print("Command errors:", e.stderr.decode())
