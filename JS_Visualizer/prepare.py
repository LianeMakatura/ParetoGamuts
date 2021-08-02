import zipfile
import argparse

json_file_paths = [
    "vis/data/bike/Rocker/fronts.json",
    "vis/data/bike/Stay/fronts.json",
    "vis/data/dense_Lamp_Stability_Mass_Focal_Point/fronts.json",
    "vis/data/keep_Turbine_Mass_Power/designExp.json",
    "vis/data/keep_Turbine_Mass_Power/fronts.json",
    "vis/data/keep100_Gridshell_nonsymb_Morning_Power_Output_Evening_Power_Output/fronts.json"
]

def compress_file(p):
    with zipfile.ZipFile(p + '.zip', 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as myzip:
        myzip.write(p)

def decompress_file(p):
    with zipfile.ZipFile(p + '.zip', 'r') as zip_file:
        zip_file.extractall("")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action")
    args = parser.parse_args()

    apply = lambda f: [f(p) for p in json_file_paths]

    if args.action == 'compress':
        apply(compress_file)
    elif args.action == 'decompress' or args.action == '':
        apply(decompress_file)

if __name__ == "__main__":
    main()