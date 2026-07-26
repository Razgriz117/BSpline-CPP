import os
import sys
import yaml
import argparse
import subprocess


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--file", type=str, default="input.yml", help="YML file with input values")

    args = parser.parse_args()
    return args


def load_input(input_path: str) -> dict:

    input = dict()

    # load input from file
    data = dict()
    with open(input_path, "r") as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)

    # quick data validation
    assert "potential" in data.keys(), f"Key 'potential' is missing from file {input_path}."
    assert isinstance(data["potential"], list) and len(data["potential"]) > 0, f"Key 'potential' in file {input_path} must map to a list with at least one element."
    for potential_piece in data["potential"]:
        assert isinstance(potential_piece, str), f"Elements of 'potential' in file {input_path} must be strings."

    # Each list element is already a JSON object string, e.g.
    # '{"domain": "[0, 20)", "function": "x"}'; join them into a single JSON
    # array string to pass as one argv argument to TISE executable, which parses
    # it back into a domain -> function map (see main.cpp:parsePiecewise).
    input["potential"] = "[" + ", ".join(data["potential"]) + "]"

    return input


if __name__ == "__main__":

    # example usage: /usr/bin/python3 input.py -f input1.yml

    # get input from file
    args = arguments()
    input = load_input(args.file)
    
    # move into TISE directory
    rootDir = os.getcwd()
    try:
        os.chdir("TISE")
    except FileNotFoundError:
        print("Error: Directory 'TISE' not found.")
        sys.exit(1)
    
    try:
        # delete any existing build dirs in TISE/
        subprocess.run(["rm", "-rf", "build "], check=True)

        # build
        subprocess.run(["cmake", "-S", ".", "-B", "build"], check=True)
        subprocess.run(["cmake", "--build", "build"], check=True)
        
        # return to original directory
        os.chdir(rootDir)
        print("Build completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        sys.exit(e.returncode)

    # call main.cpp with input-- each input must be a string
    # potential is passed as argv[1]: a single JSON-array string
    res = subprocess.run([
        "./TISE/build/H-BoundStates",
        input["potential"]
    ], capture_output=True, text=True)

    print(res.stdout)