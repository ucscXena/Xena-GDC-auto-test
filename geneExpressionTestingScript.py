import subprocess
import sys

if(len(sys.argv) != 2):
    print("Proper Usage: python3 geneExpressionTestingScript.py [Project Name] [Path to Project]")

projectName = sys.argv[2]
typeDict = {'fpkm': True,
            'fpkm_uq': True,
            'tpm': True,
            'star_counts': True}

typeToFile = { "fpkm": "star_fpkm",
                "fpkm_uq": "star_fpkm-uq",
                "tpm": "star_tpm",
                "star_counts": "star_counts"
                }
for type in typeDict:
    try:
        subprocess.run(["python3", "geneExpressionValidation.py", projectName, f"../{projectName}/Xena_Matrices/{projectName}.{typeToFile[type]}.tsv", type], check=True)
    except subprocess.CalledProcessError:
        typeDict[type] = False

print(typeDict)
