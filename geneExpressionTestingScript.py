import subprocess
import sys

if(len(sys.argv) < 2):
    print("Proper Usage: python3 geneExpressionTestingScript.py [Project Name]")

projects = sys.argv[1:]
results = []

print(projects)
for project in projects: 
    typeDict = {'star_counts': True,
                'tpm': True,
                'fpkm': True,
                'fpkm_uq': True,
                }
    typeToFile = {"star_counts": "star_counts",
                  "tpm": "star_tpm",
                  "fpkm": "star_fpkm",
                  "fpkm_uq": "star_fpkm-uq",
                }
    for type in typeDict:
        try:
            subprocess.run(["python3", "geneExpressionValidation.py", project, f"../{project}/Xena_Matrices/{project}.{typeToFile[type]}.tsv", type], check=True)
        except subprocess.CalledProcessError:
            typeDict[type] = False
    results.append([project, typeDict])

for result in results: 
    print(f'{result[0]}: {result[1]}')
